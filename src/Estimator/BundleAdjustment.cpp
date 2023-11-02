#include "BundleAdjustment.h"
#include "CostFunction.h"

using namespace std;

CBundleAdjustmentConfig::CBundleAdjustmentConfig()
{

}
size_t CBundleAdjustmentConfig::GetNumResiduals(const CModel& model, const CDatabase* const database) const
{
    // 计算所有已添加图像的观测数
	size_t numObservations = 0;
	for (size_t imageID : imageIDs)
	{
		const CImage& image = database->GetImage(imageID);
		numObservations += image.GetNumPoints3D(model.GetModelID());
	}

    // 计算所有不在imageIDs中的影像对该3D点的观测数
    auto GetNumObservationsForPoint = [this, &model](size_t point3DID) 
        {
        size_t numObservationsForPoint = 0;
        const CPoint3D& point3D = model.GetPoint3D(point3DID);
        for (const CTrackElement& trackElement : point3D.GetTrack().GetAllElements())
        {
            if (imageIDs.count(trackElement.imageID) == 0)
            {
                numObservationsForPoint++;
            }
        }
        return numObservationsForPoint;
        };
    for (size_t point3DID : variablePoint3DIDs)
    {
        numObservations += GetNumObservationsForPoint(point3DID);
    }
    for (size_t point3DID : constantPoint3DIDs)
    {
        numObservations += GetNumObservationsForPoint(point3DID);
    }
    return 2 * numObservations;
}

CBundleAdjuster::CBundleAdjuster(const COptions& options, const CBundleAdjustmentConfig& config, CDatabase* const database) :options(options), config(config), database(database)
{
    Check(database);
    options.CheckOptions();
    problem = nullptr;
}
bool CBundleAdjuster::Solve(CModel& model)
{
    Check(!problem);
    ceres::Problem::Options problemOptions;
    problemOptions.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    problem = new ceres::Problem(problemOptions);

    ceres::LossFunction* lossFunction = options.bundleAdjustmentOptions.CreateLossFunction();
    Setup(model, lossFunction);
    if (problem->NumResiduals() == 0)
    {
        delete lossFunction;
        delete problem;
        return false;
    }

    ceres::Solver::Options solverOptions = options.bundleAdjustmentOptions.ceresSolverOptions;
    solverOptions.minimizer_progress_to_stdout = true;
    solverOptions.logging_type = ceres::LoggingType::SILENT;
    const bool isHasSparse = (solverOptions.sparse_linear_algebra_library_type != ceres::NO_SPARSE);

    // 基于经验的选择
    const size_t maxNumImagesDirectDenseSolver = 50;
    const size_t maxNumImagesDirectSparseSolver = 1000;
    const size_t numImages = config.GetNumImages();
    if (numImages <= maxNumImagesDirectDenseSolver)
    {
        solverOptions.linear_solver_type = ceres::DENSE_SCHUR; // 密集求解器. 会计算和存储大量的密集矩阵, 在问题规模较小的时候往往是最快的
    }
    else if (numImages <= maxNumImagesDirectSparseSolver && isHasSparse)
    {
        solverOptions.linear_solver_type = ceres::SPARSE_SCHUR;
    }
    else // 间接稀疏(预条件共轭梯度)求解器
    {
        solverOptions.linear_solver_type = ceres::ITERATIVE_SCHUR;
        solverOptions.preconditioner_type = ceres::SCHUR_JACOBI;
    }

    if (problem->NumResiduals() < options.bundleAdjustmentOptions.minNumResidualsForMultiThreading)
    {
        solverOptions.num_threads = 1;
    }
    else
    {
        solverOptions.num_threads = GetEffectiveNumThreads();
    }

    std::string solverError;
    Check(solverOptions.IsValid(&solverError), solverError);
    ceres::Solve(solverOptions, problem, &summary);

    delete lossFunction;
    delete problem;
    return true;
}
void CBundleAdjuster::Setup(CModel& model, ceres::LossFunction* lossFunction)
{
    for (size_t imageID : config.GetAllImages())
    {
        AddImage(imageID, model, lossFunction);
    }
    for (size_t point3DID : config.GetAllVariablePoints())
    {
        AddPoint(point3DID, model, lossFunction);
    }
    for (size_t point3DID : config.GetAllConstantPoints())
    {
        AddPoint(point3DID, model, lossFunction);
    }
    ParameterizeCameras(model);
    ParameterizePoints(model);
}
void CBundleAdjuster::AddImage(size_t imageID, CModel& model, ceres::LossFunction* lossFunction)
{
    Check(database);
    CImage& image = database->GetImage(imageID);
    CCamera& camera = database->GetCamera(image.GetCameraID());
    const size_t modelID = model.GetModelID();

    // CostFunction假设是单位四元数
    CRigid3D& rigid3D = image.GetWorldToCamera(modelID);
    rigid3D.rotation.normalize();

    double* worldToCamera_Rotation = rigid3D.rotation.coeffs().data();
    double* worldToCamera_Translation = rigid3D.translation.data();
    double* cameraParams = camera.GetParamsData();

    const bool isConstantCameraPose = !options.bundleAdjustmentOptions.isRefineExtrinsics || config.IsConstantCameraPose(imageID);

    // 向平差问题中添加残差
    size_t numObservations = 0;
    const size_t numPoint2D = image.GetNumPoints2D();
    for (size_t point2DID = 0; point2DID < numPoint2D; point2DID++)
    {
        if (!image.IsPoint2DHasPoint3D(modelID, point2DID))
        {
            continue;
        }
        const CKeypoint& point2D = image.GetKeypoint(point2DID);
        const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
        CPoint3D& point3D = model.GetPoint3D(point3DID);
        CHECK(point3D.GetTrack().GetTrackLength() > 1);

        numObservations++;
        point3DNumObservations[point3DID]++;

        ceres::CostFunction* costFunction = nullptr;
        if (isConstantCameraPose)
        {
            costFunction = CReprojectionErrorConstantPoseCostFunction::Create(rigid3D, Eigen::Vector2d(point2D.pt.x, point2D.pt.y));
            problem->AddResidualBlock(costFunction, lossFunction, point3D.GetXYZ().data(), cameraParams);
        }
        else
        {
            costFunction = CReprojectionErrorCostFunction::Create(Eigen::Vector2d(point2D.pt.x, point2D.pt.y));
            problem->AddResidualBlock(costFunction, lossFunction, worldToCamera_Rotation, worldToCamera_Translation, point3D.GetXYZ().data(), cameraParams);
        }
    }

    if (numObservations > 0)
    {
        cameraIDs.insert(image.GetCameraID());
        if (!isConstantCameraPose)
        {
            SetQuaternionManifold(problem, worldToCamera_Rotation);
            if (config.IsConstantCameraPositions(imageID)) 
            {
                const vector<int>& constant_position_idxs = config.GetConstantCameraPositions(imageID);
                SetSubsetManifold(3, constant_position_idxs, problem, worldToCamera_Translation);
            }
        }
    }
}
void CBundleAdjuster::AddPoint(size_t point3DID, CModel& model, ceres::LossFunction* lossFunction)
{
    Check(database);
    CPoint3D& point3D = model.GetPoint3D(point3DID);
    const vector<CTrackElement>& trackElements = point3D.GetTrack().GetAllElements();
    const size_t modelID = model.GetModelID();

    // 3D点是否已经完全包含在问题中? 也就是说, 它的整个轨迹是否都已经被包含
    if (point3DNumObservations[point3DID] == trackElements.size())
    {
        return;
    }
    for (const CTrackElement& trackElement : trackElements)
    {
        if (config.IsExistImage(trackElement.imageID))
        {
            continue;
        }
        point3DNumObservations[point3DID]++;

        CImage& image = database->GetImage(trackElement.imageID);
        const size_t cameraID = image.GetCameraID();
        CCamera& camera = database->GetCamera(cameraID);
        const CKeypoint& point2D = image.GetKeypoint(trackElement.point2DIndex);

        image.GetWorldToCamera(modelID).rotation.normalize();
        if (cameraIDs.count(cameraID) == 0)
        {
            cameraIDs.insert(cameraID);
            config.SetConstantCameraIntrinsics(cameraID);
        }
        ceres::CostFunction* costFunction = CReprojectionErrorConstantPoseCostFunction::Create(image.GetWorldToCamera(modelID), Eigen::Vector2d(point2D.pt.x, point2D.pt.y));
        problem->AddResidualBlock(costFunction, lossFunction, point3D.GetXYZ().data(), camera.GetParamsData());
    }
}
void CBundleAdjuster::ParameterizeCameras(CModel& model)
{
    Check(database);
    const bool isConstantCamera = !options.bundleAdjustmentOptions.isRefineFocalLength && !options.bundleAdjustmentOptions.isRefinePrincipalPoint && !options.bundleAdjustmentOptions.isRefineExtraParams;
    for (size_t cameraID : cameraIDs)
    {
        CCamera& camera = database->GetCamera(cameraID);

        if (isConstantCamera || config.IsConstantCameraIntrinsics(cameraID))
        {
            problem->SetParameterBlockConstant(camera.GetParamsData());
            continue;
        }
        else
        {
            vector<int> constCameraParams;
            if (!options.bundleAdjustmentOptions.isRefineFocalLength)
            {
                size_t focalLengthIndex_X, focalLengthIndex_Y;
                camera.GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
                constCameraParams.push_back(focalLengthIndex_X);
                if (focalLengthIndex_X != focalLengthIndex_Y)
                {
                    constCameraParams.push_back(focalLengthIndex_Y);
                }
            }
            if (!options.bundleAdjustmentOptions.isRefinePrincipalPoint)
            {
                size_t principalPointIndex_X, principalPointIndex_Y;
                camera.GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
                constCameraParams.push_back(principalPointIndex_X);
                if (principalPointIndex_X != principalPointIndex_Y)
                {
                    constCameraParams.push_back(principalPointIndex_Y);
                }
            }
            if (!options.bundleAdjustmentOptions.isRefineExtraParams)
            {
                vector<size_t> paramsIndex = camera.GetExtraParamsIndex();
                constCameraParams.insert(constCameraParams.end(), paramsIndex.begin(), paramsIndex.end());
            }
            if (!constCameraParams.empty())
            {
                SetSubsetManifold(camera.GetParamsNum(), constCameraParams, problem, camera.GetParamsData());
            }
        }
    }
}
void CBundleAdjuster::ParameterizePoints(CModel& model)
{
    for (const auto& pair : point3DNumObservations)
    {
        CPoint3D& point3D = model.GetPoint3D(pair.first);
        if (point3D.GetTrack().GetTrackLength() > pair.second)
        {
            problem->SetParameterBlockConstant(point3D.GetXYZ().data());
        }
    }
    for (size_t point3DID : config.GetAllConstantPoints())
    {
        CPoint3D& point3D = model.GetPoint3D(point3DID);
        problem->SetParameterBlockConstant(point3D.GetXYZ().data());
    }
}
