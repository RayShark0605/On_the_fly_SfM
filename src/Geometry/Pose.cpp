#include "Pose.h"

#include "Triangulation.h"
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
    const double minDepth = std::numeric_limits<double>::epsilon();
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