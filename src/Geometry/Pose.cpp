#include "Pose.h"

#include "Triangulation.h"
#include "../Estimator/CostFunction.h"
using namespace std;

Eigen::Matrix3d ComputeClosestRotationMatrix(const Eigen::Matrix3d& matrix)
{
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * (svd.matrixV().transpose());
    if (R.determinant() < 0.0) 
    {
        R *= -1.0;
    }
    return R;
}
bool DecomposeProjectionMatrix(const Eigen::Matrix3x4d& projectionMatrix, Eigen::Matrix3d& K, Eigen::Matrix3d& R, Eigen::Vector3d& T)
{
    Eigen::Matrix3d RR;
    Eigen::Matrix3d QQ;
    DecomposeMatrixRQ(projectionMatrix.leftCols<3>().eval(), RR, QQ);

    R = ComputeClosestRotationMatrix(QQ);

    const double det_K = RR.determinant();
    if (det_K == 0)
    {
        return false;
    }
    else if (det_K > 0)
    {
        K = RR;
    }
    else
    {
        K = -RR;
    }

    for (int i = 0; i < 3; i++)
    {
        if (K(i, i) < 0.0)
        {
            K.col(i) = -K.col(i);
            R.row(i) = -R.row(i);
        }
    }
    T = K.triangularView<Eigen::Upper>().solve(projectionMatrix.col(3));
    if (det_K < 0)
    {
        T = -T;
    }
    return true;
}
Eigen::Matrix3d CrossProductMatrix(const Eigen::Vector3d& vector)
{
    Eigen::Matrix3d matrix;
    matrix << 0, -vector(2), vector(1), vector(2), 0, -vector(0), -vector(1), vector(0), 0;
    return matrix;
}
void RotationMatrixToEulerAngles(const Eigen::Matrix3d& R, double& rx, double& ry, double& rz)
{
    rx = atan2(R(2, 1), R(2, 2));
    ry = asin(-R(2, 0));
    rz = atan2(R(1, 0), R(0, 0));

    rx = isnan(rx) ? 0 : rx;
    ry = isnan(ry) ? 0 : ry;
    rz = isnan(rz) ? 0 : rz;
}
Eigen::Matrix3d EulerAnglesToRotationMatrix(double rx, double ry, double rz)
{
    const Eigen::Matrix3d Rx = Eigen::AngleAxisd(rx, Eigen::Vector3d::UnitX()).toRotationMatrix();
    const Eigen::Matrix3d Ry = Eigen::AngleAxisd(ry, Eigen::Vector3d::UnitY()).toRotationMatrix();
    const Eigen::Matrix3d Rz = Eigen::AngleAxisd(rz, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    return Rz * Ry * Rx;
}
Eigen::Quaterniond AverageQuaternions(const vector<Eigen::Quaterniond>& quats, const vector<double>& weights)
{
    Check(quats.size() == weights.size());
    Check(!quats.empty());

    if (quats.size() == 1) 
    {
        return quats[0];
    }

    Eigen::Matrix4d A = Eigen::Matrix4d::Zero();
    double weightSum = 0;

    for (size_t i = 0; i < quats.size(); i++) 
    {
        Check(weights[i] > 0);
        const Eigen::Vector4d qvec = quats[i].normalized().coeffs();
        A += weights[i] * qvec * qvec.transpose();
        weightSum += weights[i];
    }

    A.array() /= weightSum;

    const Eigen::Matrix4d eigenVectors = Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d>(A).eigenvectors();
    const Eigen::Vector4d averageQvec = eigenVectors.col(3);
    return Eigen::Quaterniond(averageQvec(3), averageQvec(0), averageQvec(1), averageQvec(2));
}
CRigid3D InterpolateCameraPoses(const CRigid3D& worldToCamera1, const CRigid3D& worldToCamera2, double t)
{
    const Eigen::Vector3d translation12 = worldToCamera2.translation - worldToCamera1.translation;
    return CRigid3D(worldToCamera1.rotation.slerp(t, worldToCamera2.rotation), worldToCamera1.translation + translation12 * t);
}
double CalculateDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D)
{
    const double projectedZ = worldToCamera.row(2).dot(point3D.homogeneous());
    return projectedZ * worldToCamera.col(2).norm();
}
bool CheckCheirality(const Eigen::Matrix3d& R, const Eigen::Vector3d& t, const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2, vector<Eigen::Vector3d>& points3D)
{
    Check(points1.size() == points2.size());
    const Eigen::Matrix3x4d projectionMatrix1 = Eigen::Matrix3x4d::Identity();
    Eigen::Matrix3x4d projectionMatrix2;
    projectionMatrix2.leftCols<3>() = R;
    projectionMatrix2.col(3) = t;
    const double minDepth = numeric_limits<double>::epsilon();
    const double maxDepth = 1000.0f * (R.transpose() * t).norm();
    points3D.clear();
    for (size_t i = 0; i < points1.size(); i++)
    {
        const Eigen::Vector3d point3D = TriangulatePoint(projectionMatrix1, projectionMatrix2, points1[i], points2[i]);
        const double depth1 = CalculateDepth(projectionMatrix1, point3D);
        if (depth1 > minDepth && depth1 < maxDepth)
        {
            const double depth2 = CalculateDepth(projectionMatrix2, point3D);
            if (depth2 > minDepth && depth2 < maxDepth)
            {
                points3D.push_back(point3D);
            }
        }
    }
    return !points3D.empty();
}
CRigid3D TransformCameraWorld(const CSim3D& oldWorldToNewWorld, const CRigid3D& oldWorldToCamera)
{
    const CSim3D newWorldToCamera = CSim3D(1, oldWorldToCamera.rotation, oldWorldToCamera.translation) * oldWorldToNewWorld.Inverse();
    return CRigid3D(newWorldToCamera.rotation, newWorldToCamera.translation * oldWorldToNewWorld.scale);
}
bool EstimateAbsolutePose(const COptions& options, const vector<Eigen::Vector2d>& points2D, const vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, size_t& numInliers, vector<char>& inlierMask)
{
    options.CheckOptions();

    vector<double> focalLengthFactors;
    if (options.reconstructionOptions.isEstimateAbsPoseFocalLength)
    {
        // 使用二次函数生成焦距因子, 以便对较小的焦距进行更多的样本抽取
        focalLengthFactors.reserve(options.reconstructionOptions.absPoseNumFocalLengthSamples + 1);
        const double fstep = 1.0 / options.reconstructionOptions.absPoseNumFocalLengthSamples;
        const double fscale = options.reconstructionOptions.absPoseMaxFocalLengthRatio - options.reconstructionOptions.absPoseMinFocalLengthRatio;
        double focal = 0;
        for (size_t i = 0; i <= options.reconstructionOptions.absPoseNumFocalLengthSamples; i++, focal += fstep)
        {
            focalLengthFactors.push_back(options.reconstructionOptions.absPoseMinFocalLengthRatio + fscale * focal * focal);
        }
    }
    else
    {
        focalLengthFactors.reserve(1);
        focalLengthFactors.push_back(1);
    }

    size_t focalLengthIndex_X, focalLengthIndex_Y;
    camera.GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);

    vector<CP3PEstimateRANSACReport> reports(focalLengthFactors.size());
#pragma omp parallel for
    for (int i = 0; i < focalLengthFactors.size(); i++)
    {
        CCamera scaledCamera = camera;
        scaledCamera.GetParams(focalLengthIndex_X) *= focalLengthFactors[i];
        if (focalLengthIndex_X != focalLengthIndex_Y)
        {
            scaledCamera.GetParams(focalLengthIndex_Y) *= focalLengthFactors[i];
        }

        vector<Eigen::Vector2d> points2DInCamera(points2D.size());
        for (size_t j = 0; j < points2D.size(); j++)
        {
            points2DInCamera[j] = scaledCamera.ImageToCamera(points2D[j]);
        }

        COptions customOptions = options;
        customOptions.reconstructionOptions.absPoseRANSACOoptions.maxError = scaledCamera.ImageToCameraThreshold(options.reconstructionOptions.absPoseRANSACOoptions.maxError);

        CP3PEstimator p3pEstimator;
        reports[i] = p3pEstimator.EstimateLoRANSAC(points2DInCamera, points3D, customOptions.reconstructionOptions.absPoseRANSACOoptions);
    }

    numInliers = 0;
    size_t bestModelIndex = numeric_limits<size_t>::max();
    for (size_t i = 0; i < focalLengthFactors.size(); i++)
    {
        if (reports[i].isSuccess && reports[i].support.numInliers > numInliers)
        {
            numInliers = reports[i].support.numInliers;
            bestModelIndex = i;
        }
    }

    if (numInliers == 0)
    {
        Check(bestModelIndex == numeric_limits<size_t>::max());
        return false;
    }
    inlierMask = reports[bestModelIndex].inlierMask;
    const double focalLengthFactor = focalLengthFactors[bestModelIndex];
    const Eigen::Matrix3x4d worldToCameraMatrix = reports[bestModelIndex].model;

    if (options.reconstructionOptions.isEstimateAbsPoseFocalLength)
    {
        camera.GetParams(focalLengthIndex_X) *= focalLengthFactor;
        if (focalLengthIndex_X != focalLengthIndex_Y)
        {
            camera.GetParams(focalLengthIndex_Y) *= focalLengthFactor;
        }
    }
    worldToCamera.rotation = worldToCameraMatrix.leftCols<3>();
    worldToCamera.translation = worldToCameraMatrix.col(3);
    if (worldToCamera.rotation.coeffs().array().isNaN().any() || worldToCamera.translation.array().isNaN().any())
    {
        return false;
    }
    return true;
}
bool RefineAbsolutePose(const COptions& options, const vector<char>& inlierMask, const vector<Eigen::Vector2d>& points2D, const vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, Eigen::Matrix6d* worldToCameraCov = nullptr)
{
    Check(inlierMask.size() == points2D.size() && points2D.size() == points3D.size());
    options.CheckOptions();

    ceres::CauchyLoss* lossFunction = new ceres::CauchyLoss(options.reconstructionOptions.absPoseRefineLossFunctionScale);
    double* cameraParams = camera.GetParamsData();
    double* worldToCameraRotation = worldToCamera.rotation.coeffs().data();
    double* worldToCameraTranslation = worldToCamera.translation.data();

    ceres::Problem::Options problemOptions;
    problemOptions.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problemOptions);
    for (size_t i = 0; i < points2D.size(); i++)
    {
        if (!inlierMask[i])
        {
            continue;
        }
        ceres::CostFunction* costFunction = CReprojectionErrorConstantPoint3DCostFunction::Create(points2D[i], points3D[i]);
        problem.AddResidualBlock(costFunction, lossFunction, worldToCameraRotation, worldToCameraTranslation, cameraParams);
    }
    if (problem.NumResiduals() > 0)
    {
        SetQuaternionManifold(&problem, worldToCameraRotation);
        if (!options.reconstructionOptions.isRefineAbsPoseFocalLength && !options.reconstructionOptions.isRefineAbsPoseExtraParams)
        {
            problem.SetParameterBlockConstant(camera.GetParamsData());
        }
        else
        {
            // 总是设置主点为固定的
            vector<int> cameraParamsConst;
            size_t principalPointIndexX, principalPointIndexY;
            camera.GetPrincipalPointIndex(principalPointIndexX, principalPointIndexY);
            cameraParamsConst.push_back(principalPointIndexX);
            if (principalPointIndexX != principalPointIndexY)
            {
                cameraParamsConst.push_back(principalPointIndexY);
            }

            if (!options.reconstructionOptions.isRefineAbsPoseFocalLength)
            {
                size_t focalLengthIndexX, focalLengthIndexY;
                camera.GetFocalLengthIndex(focalLengthIndexX, focalLengthIndexY);
                cameraParamsConst.push_back(focalLengthIndexX);

                if (focalLengthIndexX != focalLengthIndexY)
                {
                    cameraParamsConst.push_back(focalLengthIndexY);
                }
            }

            if (!options.reconstructionOptions.isRefineAbsPoseExtraParams)
            {
                const vector<size_t> extraParamsIndex = camera.GetExtraParamsIndex();
                cameraParamsConst.insert(cameraParamsConst.end(), extraParamsIndex.begin(), extraParamsIndex.end());
            }

            if (cameraParamsConst.size() == camera.GetParamsNum())
            {
                problem.SetParameterBlockConstant(camera.GetParamsData());
            }
            else
            {
                SetSubsetManifold(camera.GetParamsNum(), cameraParamsConst, &problem, camera.GetParamsData());
            }
        }
    }

    ceres::Solver::Options solverOptions;
    solverOptions.gradient_tolerance = options.reconstructionOptions.absPoseRefineGradientTolerance;
    solverOptions.max_num_iterations = options.reconstructionOptions.absPoseRefineMaxNumIterations;
    solverOptions.linear_solver_type = ceres::DENSE_QR;
    solverOptions.num_threads = 1;
    solverOptions.logging_type = ceres::SILENT;
    solverOptions.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(solverOptions, &problem, &summary);

    if (problem.NumResiduals() > 0 && worldToCameraCov)
    {
        ceres::Covariance::Options options;
        ceres::Covariance covariance(options);
        vector<const double*> parameterBlocks = { worldToCameraRotation, worldToCameraTranslation };
        if (!covariance.Compute(parameterBlocks, &problem)) 
        {
            return false;
        }
        // 旋转协方差是在四元数的切线空间中估算的, 这与 3-DoF 轴角局部参数化相对应
        covariance.GetCovarianceMatrixInTangentSpace(parameterBlocks, worldToCameraCov->data());
    }
    return summary.IsSolutionUsable();
}








