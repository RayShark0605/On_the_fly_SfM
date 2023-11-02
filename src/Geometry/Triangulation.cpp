#include "Triangulation.h"
#include "../Scene/Database.h"
#include "../Scene/Model.h"
#include "../Scene/Image.h"

using namespace std;

Eigen::Vector3d TriangulatePoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
{
	// Step 1. 根据双视几何之间的三角测量方程来构建矩阵A
	Eigen::Matrix4d A;
	A.row(0) = point1(0) * worldToCamera1.row(2) - worldToCamera1.row(0);
	A.row(1) = point1(1) * worldToCamera1.row(2) - worldToCamera1.row(1);
	A.row(2) = point2(0) * worldToCamera2.row(2) - worldToCamera2.row(0);
	A.row(3) = point2(1) * worldToCamera2.row(2) - worldToCamera2.row(1);

	// Step 2. 对A做奇异值分解, 计算其所有右奇异向量
	Eigen::JacobiSVD<Eigen::Matrix4d> svd(A, Eigen::ComputeFullV);

	// Step 3. 最后一个右奇异向量对应于矩阵A的最小奇异值, 在三角测量问题中, 这个向量就表示三维点. 通过hnormalized转为齐次坐标
	return svd.matrixV().col(3).hnormalized();
}
Eigen::Vector3d TriangulateOptimalPoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
{
	const CRigid3D worldToCamera1_Rigid3D(Eigen::Quaterniond(worldToCamera1.leftCols<3>()), worldToCamera1.col(3));
	const CRigid3D worldToCamera2_Rigid3D(Eigen::Quaterniond(worldToCamera2.leftCols<3>()), worldToCamera2.col(3));
	const CRigid3D camera1ToCamera2_Rigid3D = worldToCamera2_Rigid3D * worldToCamera1_Rigid3D.Inverse();
	const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2_Rigid3D);

	Eigen::Vector2d optimalPoint1;
	Eigen::Vector2d optimalPoint2;
	FindOptimalImageObservations(E, point1, point2, optimalPoint1, optimalPoint2);
	return TriangulatePoint(worldToCamera1, worldToCamera2, optimalPoint1, optimalPoint2);
}
Eigen::Vector3d TriangulateMultiViewPoint(const vector<Eigen::Matrix3x4d>& worldToCameras, const vector<Eigen::Vector2d>& points2D)
{
	Check(worldToCameras.size() == points2D.size());

	// Step 1. 初始化4×4的矩阵A, 它将用于累积所有视图的信息
	Eigen::Matrix4d A = Eigen::Matrix4d::Zero();

	// Step 2. 遍历所有的视图和每个视图中对应的2D点
	for (size_t i = 0; i < points2D.size(); i++)
	{
		// Step 2.1. 将2D点转换为齐次坐标, 并进行归一化
		const Eigen::Vector3d point = points2D[i].homogeneous().normalized();

		// Step 2.2. 计算误差项(误差矩阵): 原始投影矩阵减去调整后的投影矩阵
		const Eigen::Matrix3x4d term = worldToCameras[i] - point * point.transpose() * worldToCameras[i];

		// Step 2.3. 对误差项做平方求和
		A += term.transpose() * term;
	}

	// Step 3. 使用特征值分解来求最小二乘问题. 这里使用的是自伴矩阵(self-adjoint)的特征值求解器, 它数值更稳定
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigenSolver(A);

	// Step 4. 返回与最小特征值对应的特征向量, 它就是最优的3D点位置
	return eigenSolver.eigenvectors().col(0).hnormalized();
}
vector<Eigen::Vector3d> TriangulatePoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2)
{
	Check(points1.size() == points2.size());
	vector<Eigen::Vector3d> points3D(points1.size());
	for (size_t i = 0; i < points3D.size(); i++) 
	{
		points3D[i] = TriangulatePoint(worldToCamera1, worldToCamera2, points1[i], points2[i]);
	}
	return points3D;
}
vector<Eigen::Vector3d> TriangulateOptimalPoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2)
{
	Check(points1.size() == points2.size());
	vector<Eigen::Vector3d> points3D(points1.size());
	for (size_t i = 0; i < points3D.size(); i++)
	{
		points3D[i] = TriangulateOptimalPoint(worldToCamera1, worldToCamera2, points1[i], points2[i]);
	}
	return points3D;
}
double CalculateTriangulationAngle(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const Eigen::Vector3d& point3D)
{
	const double baselineLengthSquared = (projectionCenter1 - projectionCenter2).squaredNorm();

	const double rayLengthSquared1 = (point3D - projectionCenter1).squaredNorm();
	const double rayLengthSquared2 = (point3D - projectionCenter2).squaredNorm();

	// 使用余弦定理计算交会角
	const double denominator = 2.0 * sqrt(rayLengthSquared1 * rayLengthSquared2);
	if (denominator == 0.0)
	{
		return 0.0;
	}
	const double nominator = rayLengthSquared1 + rayLengthSquared2 - baselineLengthSquared;
	const double angle = abs(acos(nominator / denominator));

	// 三角测量在锐角(比较远的同名点)和钝角(比较近的同名点)的情况下是不稳定的, 所以总是要计算两条相交射线形成的两种角度(锐角和钝角)中的较小值
	return min(angle, M_PI - angle);
}
vector<double> CalculateTriangulationAngles(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const vector<Eigen::Vector3d>& points3D)
{
	const double baselineSquared = (projectionCenter1 - projectionCenter2).squaredNorm();
	vector<double> angles(points3D.size());
	for (size_t i = 0; i < points3D.size(); i++)
	{
		const double rayLengthSquared1 = (points3D[i] - projectionCenter1).squaredNorm();
		const double rayLengthSquared2 = (points3D[i] - projectionCenter2).squaredNorm();

		const double denominator = 2.0 * sqrt(rayLengthSquared1 * rayLengthSquared2);
		if (denominator == 0.0)
		{
			angles[i] = 0.0;
			continue;
		}
		const double nominator = rayLengthSquared1 + rayLengthSquared2 - baselineSquared;
		const double angle = abs(acos(nominator / denominator));
		angles[i] = min(angle, M_PI - angle);
	}
	return angles;
}
bool EstimateTriangulation(const CTriangulationOptions& options, const vector<CTriangulationPoint>& points, const vector<CTriangulationPose>& poses, vector<char>& inlierMask, Eigen::Vector3d& XYZ)
{
	Check(points.size() >= 2 && points.size() == poses.size());
	options.CheckOptions();

	CTriangulationEstimator triangulationEstimator(options.minTriAngle_Deg * M_PI / 180, options.residualType);
	CTriangulationEstimator triangulationEstimator_Local(options.minTriAngle_Deg * M_PI / 180, options.residualType);
	CInlierSupportMeasurer inlierSupportMeasurer;
	CCombinationSampler combinationSampler(triangulationEstimator.minNumSamples);
	CLORANSAC loransac(options.ransacOptions, &triangulationEstimator, &triangulationEstimator_Local, &inlierSupportMeasurer, &combinationSampler);
	const CRANSACReport report = loransac.Estimate<CTriangulationPoint, CTriangulationPose, Eigen::Vector3d>(points, poses);
	if (!report.isSuccess) 
	{
		return false;
	}
	inlierMask = report.inlierMask;
	XYZ = report.model;
	return true;
}


CTriangulator::CTriangulator(CDatabase* const database, CModel& model) : database(database), model(model), modelID(model.GetModelID())
{
	Check(database);
}
size_t CTriangulator::TriangulateImage(const COptions& options, size_t imageID)
{
	options.CheckOptions();
	ClearCaches();

	const CImage& image = database->GetImage(imageID);
	bool isRegistered = image.IsRegistered(modelID);
	Check(model.IsImageRegistered(imageID) == isRegistered);
	if (!isRegistered)
	{
		return 0;
	}

	const CCamera& camera = database->GetCamera(image.GetCameraID());
	bool isBogusParams = IsCameraBogusParams(options, camera);
	Check(camera.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam) == isBogusParams);
	if (isBogusParams)
	{
		return 0;
	}

	CCorrData refCorrData;
	refCorrData.imageID = imageID;
	refCorrData.image = &image;
	refCorrData.camera = &camera;

	vector<CCorrData> corrsData;
	size_t numTris = 0;
	for (size_t point2DID = 0; point2DID < image.GetNumPoints2D(); point2DID++)
	{
		const size_t numTriangulated = Find(options, imageID, point2DID, options.estimateTriangulationOptions.maxTransitivity, corrsData);
		if (corrsData.empty())
		{
			continue;
		}
		
		const CKeypoint& point2D = image.GetKeypoint(point2DID);
		refCorrData.point2DID = point2DID;
		refCorrData.point2D = &point2D;

		if (numTriangulated != 0)
		{
			numTris += Continue(options, refCorrData, corrsData);
		}
		corrsData.push_back(refCorrData);
		numTris += Create(options, corrsData);
	}
	return numTris;
}
size_t CTriangulator::CompleteImage(const COptions& options, size_t imageID)
{
	options.CheckOptions();
	ClearCaches();

	const CImage& image = database->GetImage(imageID);
	bool isRegistered = image.IsRegistered(modelID);
	Check(model.IsImageRegistered(imageID) == isRegistered);
	if (!isRegistered)
	{
		return 0;
	}

	const CCamera& camera = database->GetCamera(image.GetCameraID());
	bool isBogusParams = IsCameraBogusParams(options, camera);
	Check(camera.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam) == isBogusParams);
	if (isBogusParams)
	{
		return 0;
	}

	CTriangulationOptions triOptions;
	triOptions.residualType = CTriangulationResidualType::CReprojectionError;
	triOptions.ransacOptions.maxError = options.estimateTriangulationOptions.completeMaxReprojectionError;
	triOptions.ransacOptions.confidence = 0.9999;
	triOptions.ransacOptions.minInlierRatio = 0.02;
	triOptions.ransacOptions.maxNumTrials = 10000;

	CCorrData refCorrData;
	refCorrData.imageID = imageID;
	refCorrData.image = &image;
	refCorrData.camera = &camera;

	vector<CCorrData> corrsData;
	size_t numTris = 0;
	for (size_t point2DID = 0; point2DID < image.GetNumPoints2D(); point2DID++)
	{
		if (image.IsPoint2DHasPoint3D(modelID, point2DID))
		{
			numTris += Complete(options, image.GetPoint3DID(point2DID, modelID));
			continue;
		}
		if (options.estimateTriangulationOptions.isIgnoreTwoViewTracks && IsTwoViewObservation(image, point2DID))
		{
			continue;
		}
		const size_t numTriangulated = Find(options, imageID, point2DID, options.estimateTriangulationOptions.maxTransitivity, corrsData);
		if (numTriangulated || corrsData.empty())
		{
			continue;
		}
		const CKeypoint& point2D = image.GetKeypoint(point2DID);
		refCorrData.point2DID = point2DID;
		refCorrData.point2D = &point2D;
		corrsData.push_back(refCorrData);

		vector<CTriangulationPoint> pointData;
		vector<CTriangulationPose> poseData;
		pointData.resize(corrsData.size());
		poseData.resize(corrsData.size());
		for (size_t i = 0; i < corrsData.size(); i++)
		{
			const CCorrData& corrData = corrsData[i];
			pointData[i].point.x() = corrData.point2D->pt.x;
			pointData[i].point.y() = corrData.point2D->pt.y;
			pointData[i].pointNormalized = corrData.camera->ImageToCamera(pointData[i].point);

			poseData[i].projectionCenter = corrData.image->GetProjectionCenter(modelID);
			poseData[i].projectionMatrix = corrData.image->GetWorldToCamera(modelID).ToMatrix();
			poseData[i].camera = corrData.camera;
		}
		if (pointData.size() <= 15)
		{
			triOptions.ransacOptions.minNumTrials = NChooseK(pointData.size(), 2);
		}
		Eigen::Vector3d XYZ;
		vector<char> inlierMask;
		if (!EstimateTriangulation(triOptions, pointData, poseData, inlierMask, XYZ))
		{
			continue;
		}
		
		CTrack track;
		track.Reserve(corrsData.size());
		for (size_t i = 0; i < inlierMask.size(); i++)
		{
			if (inlierMask[i])
			{
				const CCorrData& corrData = corrsData[i];
				track.AddElement(corrData.imageID, corrData.point2DID);
				numTris++;
			}
		}

		const size_t point3DID = model.AddPoint3D(XYZ, move(track));
		modifiedPoint3DIDs.insert(point3DID);
	}
	return numTris;
}
size_t CTriangulator::CompleteTracks(const COptions& options, const unordered_set<size_t>& points3DIDs)
{
	options.CheckOptions();
	ClearCaches();

	size_t numCompleted = 0;
	for (size_t point3DID : points3DIDs)
	{
		numCompleted += Complete(options, point3DID);
	}
	return numCompleted;
}
size_t CTriangulator::CompleteAllTracks(const COptions& options)
{
	options.CheckOptions();
	ClearCaches();

	const unordered_set<size_t> points3DIDs = model.GetAllPoints3D();
	size_t numCompleted = 0;
	for (size_t point3DID : points3DIDs)
	{
		numCompleted += Complete(options, point3DID);
	}
	return numCompleted;
}
size_t CTriangulator::MergeTracks(const COptions& options, const unordered_set<size_t>& points3DIDs)
{
	options.CheckOptions();
	ClearCaches();

	size_t numMerged = 0;
	for (size_t point3DID : points3DIDs)
	{
		numMerged += Merge(options, point3DID);
	}
	return numMerged;
}
size_t CTriangulator::MergeAllTracks(const COptions& options)
{
	options.CheckOptions();
	ClearCaches();

	const unordered_set<size_t> points3DIDs = model.GetAllPoints3D();
	size_t numMerged = 0;
	for (size_t point3DID : points3DIDs)
	{
		numMerged += Merge(options, point3DID);
	}
	return numMerged;
}
size_t CTriangulator::Retriangulate(const COptions& options)
{
	options.CheckOptions();
	ClearCaches();

	COptions retriangulateOptions = options;
	retriangulateOptions.estimateTriangulationOptions.continueMaxAngleError = options.estimateTriangulationOptions.retriangulateMaxAngleError;

	size_t numTris = 0;
	const CImagePairsType& imagePairs = model.GetAllImagePairs();
	for (const auto& pair : imagePairs)
	{
		const std::pair<size_t, size_t>& imagePair = pair.first;
		const CImagePairStatus& imagePairStatus = pair.second;
		Check(imagePair.first < imagePair.second);
		const double triRatio = imagePairStatus.numTriCorrs * 1.0 / imagePairStatus.numTotalCorrs;
		if (triRatio >= options.estimateTriangulationOptions.retriangulateMinRatio)
		{
			continue;
		}
		if (!model.IsImageRegistered(imagePair.first) || !model.IsImageRegistered(imagePair.second))
		{
			continue;
		}
		if (retriangulateNumTrials[imagePair] >= options.estimateTriangulationOptions.retriangulateMaxTrials)
		{
			continue;
		}
		retriangulateNumTrials[imagePair]++;

		const CImage image1 = database->GetImage(imagePair.first);
		const CImage image2 = database->GetImage(imagePair.second);
		Check(image1.IsRegistered(modelID) && image2.IsRegistered(modelID));
		const CCamera camera1 = database->GetCamera(image1.GetCameraID());
		const CCamera camera2 = database->GetCamera(image2.GetCameraID());
		if (IsCameraBogusParams(options, camera1) || IsCameraBogusParams(options, camera2))
		{
			continue;
		}
		Check(!camera1.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam));
		Check(!camera2.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam));

		const CSIFTMatches& corrs = database->GetTwoViewGeometry(imagePair.first, imagePair.second).inlierMatches;
		for (const CSIFTMatch& match : corrs)
		{
			// 这里有两种可能的情况: 两个点要么属于同一个3D点, 要么属于不同的3D点. 
			// 在前一种情况下, 没有什么需要做的. 在后一种情况下, 我们不尝试重新前方交会, 因为重新前方交会出的对应关系很可能是不准确的, 因此如果合并, 将破坏两个3D点
			bool isHasPoint3D1 = image1.IsPoint2DHasPoint3D(modelID, match.point2DIndex1);
			bool isHasPoint3D2 = image2.IsPoint2DHasPoint3D(modelID, match.point2DIndex2);
			if (isHasPoint3D1 && isHasPoint3D2)
			{
				continue;
			}

			const CKeypoint& point2D1 = image1.GetKeypoint(match.point2DIndex1);
			const CKeypoint& point2D2 = image2.GetKeypoint(match.point2DIndex2);

			CCorrData corr_data1;
			corr_data1.imageID = imagePair.first;
			corr_data1.point2DID = match.point2DIndex1;
			corr_data1.image = &image1;
			corr_data1.camera = &camera1;
			corr_data1.point2D = &point2D1;

			CCorrData corr_data2;
			corr_data2.imageID = imagePair.second;
			corr_data2.point2DID = match.point2DIndex2;
			corr_data2.image = &image2;
			corr_data2.camera = &camera2;
			corr_data2.point2D = &point2D2;

			if (isHasPoint3D1 && !isHasPoint3D2)
			{
				const vector<CCorrData> corrsData = { corr_data1 };
				numTris += Continue(retriangulateOptions, corr_data2, corrsData);
			}
			else if (!isHasPoint3D1 && isHasPoint3D2)
			{
				const vector<CCorrData> corrsData = { corr_data2 };
				numTris += Continue(retriangulateOptions, corr_data1, corrsData);
			}
			else // 两个都没有3D点
			{
				const vector<CCorrData> corrsData = { corr_data1,  corr_data2 };
				numTris += Create(options, corrsData); // 不要使用较大的前方交会阈值, 因为这会在创建点时导致明显的漂移, 所以这里使用options而不是retriangulateOptions
			}
		}
	}
	return numTris;
}
void CTriangulator::AddModifiedPoint3D(size_t point3DID)
{
	modifiedPoint3DIDs.insert(point3DID);
}
const unordered_set<size_t>& CTriangulator::GetModifiedPoints3D()
{
	for (auto it = modifiedPoint3DIDs.begin(); it != modifiedPoint3DIDs.end();)
	{
		if (model.IsExistsPoint3D(*it))
		{
			it++;
		}
		else
		{
			it = modifiedPoint3DIDs.erase(it);
		}
	}
	return modifiedPoint3DIDs;
}
void CTriangulator::ClearModifiedPoints3D()
{
	modifiedPoint3DIDs.clear();
}
void CTriangulator::ClearCaches()
{
	isCameraBogusParams.clear();
	mergeTrials.clear();
}
size_t CTriangulator::Find(const COptions& options, size_t imageID, size_t point2DID, size_t transitivity, vector<CCorrData>& corrsData)
{
	CConjugatePoints corrs;
	FindTransitiveCorrespondences(imageID, point2DID, transitivity, corrs);

	corrsData.clear();
	corrsData.reserve(corrs.size());
	size_t numTriangulated = 0;
	for (const auto& pair : corrs)
	{
		if (!model.IsImageRegistered(pair.first))
		{
			continue;
		}
		const CImage& corrImage = database->GetImage(pair.first);
		const CCamera& corrCamera = database->GetCamera(corrImage.GetCameraID());
		if (IsCameraBogusParams(options, corrCamera))
		{
			continue;
		}

		CCorrData corrData;
		corrData.imageID = pair.first;
		corrData.point2DID = pair.second;
		corrData.image = &corrImage;
		corrData.camera = &corrCamera;
		corrData.point2D = &corrImage.GetKeypoint(pair.second);
		corrsData.push_back(corrData);

		if (corrImage.IsPoint2DHasPoint3D(modelID, pair.second))
		{
			numTriangulated++;
		}
	}
	return numTriangulated;
}
size_t CTriangulator::Create(const COptions& options, const vector<CCorrData>& corrsData)
{
	vector<CCorrData> createCorrsData;
	createCorrsData.reserve(corrsData.size());
	for (const CCorrData& corrData : corrsData)
	{
		if (!corrData.image->IsPoint2DHasPoint3D(modelID, corrData.point2DID))
		{
			createCorrsData.push_back(corrData);
		}
	}

	if (createCorrsData.size() < 2)
	{
		return 0;
	}
	if (createCorrsData.size() == 2 && options.estimateTriangulationOptions.isIgnoreTwoViewTracks)
	{
		const CCorrData& corrData1 = createCorrsData[0];
		if (IsTwoViewObservation(*corrData1.image, corrData1.point2DID))
		{
			return 0;
		}
	}

	vector<CTriangulationPoint> pointData(createCorrsData.size());
	vector<CTriangulationPose> poseData(createCorrsData.size());
	for (size_t i = 0; i < createCorrsData.size(); i++)
	{
		const CCorrData& corrData = createCorrsData[i];
		pointData[i].point.x() = corrData.point2D->pt.x;
		pointData[i].point.y() = corrData.point2D->pt.y;
		pointData[i].pointNormalized = corrData.camera->ImageToCamera(pointData[i].point);

		poseData[i].projectionCenter = corrData.image->GetProjectionCenter(modelID);
		poseData[i].projectionMatrix = corrData.image->GetWorldToCamera(modelID).ToMatrix();
		poseData[i].camera = corrData.camera;
	}

	CTriangulationOptions triangleOptions;
	triangleOptions.minTriAngle_Deg = options.estimateTriangulationOptions.minTriAngle_Deg;
	triangleOptions.residualType = CTriangulationResidualType::CAngularError;
	triangleOptions.ransacOptions.maxError = DegToRad(options.estimateTriangulationOptions.createMaxAngleError);
	triangleOptions.ransacOptions.confidence = 0.9999;
	triangleOptions.ransacOptions.minInlierRatio = 0.02;
	triangleOptions.ransacOptions.maxNumTrials = 10000;
	if (pointData.size() <= 15)
	{
		triangleOptions.ransacOptions.minNumTrials = NChooseK(pointData.size(), 2);
	}

	Eigen::Vector3d XYZ;
	vector<char> inlierMask;
	if (!EstimateTriangulation(triangleOptions, pointData, poseData, inlierMask, XYZ))
	{
		return 0;
	}

	CTrack track;
	track.Reserve(createCorrsData.size());
	for (size_t i = 0; i < inlierMask.size(); i++)
	{
		if (inlierMask[i])
		{
			const CCorrData& corrData = createCorrsData[i];
			track.AddElement(corrData.imageID, corrData.point2DID);
		}
	}

	const size_t trackLength = track.GetTrackLength();
	const size_t point3DID = model.AddPoint3D(XYZ, move(track));
	modifiedPoint3DIDs.insert(point3DID);

	if (createCorrsData.size() >= trackLength + 3)
	{
		return trackLength + Create(options, createCorrsData);
	}
	return trackLength;
}
size_t CTriangulator::Continue(const COptions& options, const CCorrData& refCorrData, const vector<CCorrData>& corrsData)
{
	if (refCorrData.image->IsPoint2DHasPoint3D(modelID, refCorrData.point2DID))
	{
		return 0;
	}

	double bestAngleError = numeric_limits<double>::max();
	size_t bestIndex = numeric_limits<size_t>::max();
	for (size_t index = 0; index < corrsData.size(); index++)
	{
		const CCorrData& corrData = corrsData[index];
		if (!corrData.image->IsPoint2DHasPoint3D(modelID, corrData.point2DID))
		{
			continue;
		}
		const CPoint3D& point3D = model.GetPoint3D(corrData.image->GetPoint3DID(corrData.point2DID, modelID));
		const double angleError = CalculateAngularError(Eigen::Vector2d(refCorrData.point2D->pt.x, refCorrData.point2D->pt.y), point3D.GetXYZ(), refCorrData.image->GetWorldToCamera(modelID), *refCorrData.camera);
		if (angleError < bestAngleError)
		{
			bestAngleError = angleError;
			bestIndex = index;
		}
	}

	const double maxAngleError = DegToRad(options.estimateTriangulationOptions.continueMaxAngleError);
	if (bestAngleError <= maxAngleError && bestIndex != numeric_limits<size_t>::max())
	{
		const CCorrData& corrData = corrsData[bestIndex];
		const CTrackElement trackElement(refCorrData.imageID, refCorrData.point2DID);
		Check(corrData.image->IsPoint2DHasPoint3D(modelID, corrData.point2DID));
		const size_t point3DID = corrData.image->GetPoint3DID(corrData.point2DID, modelID);
		model.AddObservation(point3DID, trackElement);
		modifiedPoint3DIDs.insert(point3DID);
		return 1;
	}
	return 0;
}
size_t CTriangulator::Merge(const COptions& options, size_t point3DID)
{
	if (!model.IsExistsPoint3D(point3DID))
	{
		return 0;
	}
	const double maxSquaredReprojectionError = options.estimateTriangulationOptions.mergeMaxReprojectionError * options.estimateTriangulationOptions.mergeMaxReprojectionError;
	const CPoint3D& point3D = model.GetPoint3D(point3DID);
	const vector<CTrackElement> trackElements = point3D.GetTrack().GetAllElements();
	for (const CTrackElement& trackElement : trackElements)
	{
		const CConjugatePoints& correspondences = database->GetImage(trackElement.imageID).GetCorrespondences(trackElement.point2DIndex).first;
		for (const auto& pair : correspondences)
		{
			if (!model.IsImageRegistered(pair.first))
			{
				continue;
			}
			const CImage& image = database->GetImage(pair.first);
			const CKeypoint& corrPoint2D = image.GetKeypoint(pair.second);
			if (!image.IsPoint2DHasPoint3D(modelID, pair.second))
			{
				continue;
			}
			const size_t corrPoint3DID = image.GetPoint3DID(pair.second, modelID);
			if (corrPoint3DID == point3DID || mergeTrials[point3DID].count(corrPoint3DID) > 0)
			{
				continue;
			}
			const CPoint3D& corrPoint3D = model.GetPoint3D(corrPoint3DID);
			mergeTrials[point3DID].insert(corrPoint3DID);
			mergeTrials[corrPoint3DID].insert(point3DID);

			const size_t point3DTrackLength = point3D.GetTrack().GetTrackLength();
			const size_t corrPoint3DTrackLength = corrPoint3D.GetTrack().GetTrackLength();

			const Eigen::Vector3d mergedXYZ = (point3DTrackLength * point3D.GetXYZ() + corrPoint3DTrackLength * corrPoint3D.GetXYZ()) / (point3DTrackLength + corrPoint3DTrackLength);
			bool isMergeSuccess = true;
			for (const CTrack* track : { &point3D.GetTrack(), &corrPoint3D.GetTrack() })
			{
				for (const CTrackElement& testTrackElement : track->GetAllElements())
				{
					const CImage& testImage = database->GetImage(testTrackElement.imageID);
					const CCamera& testCamera = database->GetCamera(testImage.GetCameraID());
					const CKeypoint& testPoint2D = testImage.GetKeypoint(testTrackElement.point2DIndex);

					if (CalculateSquaredReprojectionError(Eigen::Vector2d(testPoint2D.pt.x, testPoint2D.pt.y), mergedXYZ, testImage.GetWorldToCamera(modelID), testCamera) > maxSquaredReprojectionError)
					{
						isMergeSuccess = false;
						break;
					}
				}
				if (!isMergeSuccess)
				{
					break;
				}
			}
			if (isMergeSuccess)
			{
				const size_t numMerged = point3DTrackLength + corrPoint3DTrackLength;
				const size_t mergedPoint3DID = model.MergePoints3D(point3DID, corrPoint3DID);
				modifiedPoint3DIDs.erase(point3DID);
				modifiedPoint3DIDs.erase(corrPoint3DID);
				modifiedPoint3DIDs.insert(mergedPoint3DID);

				const size_t numMergedRecursive = Merge(options, mergedPoint3DID);
				if (numMergedRecursive > 0)
				{
					return numMergedRecursive;
				}
				return numMerged;
			}
		}
	}
	return 0;
}
size_t CTriangulator::Complete(const COptions& options, size_t point3DID)
{
	if (!model.IsExistsPoint3D(point3DID))
	{
		return 0;
	}

	size_t numCompleted = 0;
	const double maxSquaredReprojectionError = options.estimateTriangulationOptions.completeMaxReprojectionError * options.estimateTriangulationOptions.completeMaxReprojectionError;

	const CPoint3D& point3D = model.GetPoint3D(point3DID);
	vector<CTrackElement> trackElements = point3D.GetTrack().GetAllElements();

	const int maxTransitivity = options.estimateTriangulationOptions.maxTransitivity;
	for (int transitivity = 0; transitivity < maxTransitivity; transitivity++)
	{
		if (trackElements.empty())
		{
			break;
		}
		const vector<CTrackElement> preTrackElement = trackElements;
		trackElements.clear();
		for (const CTrackElement& trackElement : preTrackElement)
		{
			const CConjugatePoints& correspondences = database->GetImage(trackElement.imageID).GetCorrespondences(trackElement.point2DIndex).first;
			for (const auto& pair : correspondences)
			{
				if (!model.IsImageRegistered(pair.first))
				{
					continue;
				}
				const CImage& image = database->GetImage(pair.first);
				if (image.IsPoint2DHasPoint3D(modelID, pair.second))
				{
					continue;
				}
				const CKeypoint& point2D = image.GetKeypoint(pair.second);
				const CPoint3D& point3D = model.GetPoint3D(image.GetPoint3DID(pair.second, modelID));
				const CCamera& camera = database->GetCamera(image.GetCameraID());
				if (IsCameraBogusParams(options, camera))
				{
					continue;
				}

				if (CalculateSquaredReprojectionError(Eigen::Vector2d(point2D.pt.x, point2D.pt.y), point3D.GetXYZ(), image.GetWorldToCamera(modelID), camera) > maxSquaredReprojectionError)
				{
					continue;
				}

				const CTrackElement trackElement(pair.first, pair.second);
				model.AddObservation(point3DID, trackElement);
				modifiedPoint3DIDs.insert(point3DID);

				if (transitivity < maxTransitivity - 1)
				{
					trackElements.emplace_back(pair.first, pair.second);
				}
				numCompleted++;
			}
		}
	}
	return numCompleted;
}
bool CTriangulator::IsCameraBogusParams(const COptions& options, const CCamera& camera)
{
	const auto it = isCameraBogusParams.find(camera.GetCameraID());
	if (it != isCameraBogusParams.end())
	{
		return it->second;
	}
	const bool isBogusParams = camera.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam);
	isCameraBogusParams[camera.GetCameraID()] = isBogusParams;
	return isBogusParams;
}
bool CTriangulator::IsTwoViewObservation(const CImage& image, size_t point2DID) const
{
	const CConjugatePoints& conjugatePoints = image.GetCorrespondences(point2DID).first;
	if (conjugatePoints.size() != 1)
	{
		return false;
	}
	
	const size_t corrImageID = conjugatePoints.begin()->first;
	const size_t corrPoint2DID = conjugatePoints.begin()->second;
	const CImage& corrImage = database->GetImage(corrImageID);
	const CConjugatePoints corrConjugatePoints = corrImage.GetCorrespondences(corrPoint2DID).first;

	return (corrConjugatePoints.size() == 1);
}
void CTriangulator::FindTransitiveCorrespondences(size_t imageID, size_t point2DID, size_t transitivity, CConjugatePoints& corrs)
{
	const CImage& image = database->GetImage(imageID);
	corrs = image.GetCorrespondences(point2DID).first;

	queue<pair<size_t, size_t>> queue;
	for (const auto& pair : corrs)
	{
		queue.push({ pair.first,pair.second });
	}

	while (!queue.empty() && transitivity > 0)
	{
		const size_t currentLayerSize = queue.size();
		for (size_t i = 0; i < currentLayerSize; i++)
		{
			const size_t curImageID = queue.front().first;
			const size_t curPoint2DID = queue.front().second;
			queue.pop();
			if (corrs.find(curImageID) == corrs.end())
			{
				corrs[curImageID] = curPoint2DID;
			}
			else
			{
				Check(corrs[curImageID] == curPoint2DID);
			}
			const CConjugatePoints& moreCorrs = database->GetImage(curImageID).GetCorrespondences(curPoint2DID).first;
			for (const auto& pair : moreCorrs)
			{
				queue.push({ pair.first,pair.second });
			}
		}
		transitivity--;
	}
}
