#include "Mapper.h"
#include "BundleAdjustment.h"
#include "../Base/Math.h"

using namespace std;

void SortAndAppendNextImages(vector<pair<size_t, float>>& imageRanks, vector<size_t>& sortedImageIDs)
{
	sort(imageRanks.begin(), imageRanks.end(), [](const pair<size_t, float>& image1, const pair<size_t, float>& image2)
		{
			return image1.second > image2.second;
		});
	sortedImageIDs.reserve(sortedImageIDs.size() + imageRanks.size());
	for (const pair<size_t, float>& pair : imageRanks)
	{
		sortedImageIDs.push_back(pair.first);
	}
}

CMapper::CMapper(CModel& model, CDatabase* const database) :model(model), database(database), triangulator(new CTriangulator(database, model)), modelID(model.GetModelID())
{
	preInitImagePair = { numeric_limits<size_t>::max(), numeric_limits<size_t>::max() };
	initImagePairs.clear();
	numRegImagesPerCamera.clear();
	filteredImages.clear();
	existingImages.clear();
}
CMapper::~CMapper()
{
	if (triangulator)
	{
		delete triangulator;
	}
	initImagePairs.clear();
	numRegImagesPerCamera.clear();
	filteredImages.clear();
	existingImages.clear();
}
void CMapper::BeginReconstruction(const COptions& options)
{
	const unordered_set<size_t> regImages = model.GetAllRegImages();
	for (size_t imageID : regImages)
	{
		RegisterImageEvent(imageID);
	}
	existingImages = regImages;

	if (options.reconstructionOptions.nextImageSelectionMethod == CNextImageSelectionMethod::CMinUncertainty)
	{
		Check(database);
		const size_t numImages = database->GetImagesNum();
		for (size_t imageID = 0; imageID < numImages; imageID++)
		{
			CImage& image = database->GetImage(imageID);
			const CCamera& camera = database->GetCamera(image.GetCameraID());
			image.Setup(camera);
		}
	}
}
void CMapper::EndReconstruction(bool isDiscard)
{
	if (isDiscard)
	{
		Check(database);
		const unordered_set<size_t> regImages = model.GetAllRegImages();
		for (size_t imageID : regImages)
		{
			DeRegisterImageEvent(imageID);
			database->GetImage(imageID).DeleteModelInfo(modelID);
		}
	}
}
bool CMapper::FindInitialImagePair(const COptions& options, size_t& imageID1, size_t& imageID2)
{
	options.CheckOptions();
	imageID1 = imageID2 = numeric_limits<size_t>::max();
	const vector<size_t> findImageIDs1 = FindFirstInitialImage();
	for (size_t findImageID1 : findImageIDs1)
	{
		const vector<size_t> findImageIDs2 = FindSecondInitialImage(options, findImageID1);
		for (size_t findImageID2 : findImageIDs2)
		{
			const pair<size_t, size_t> imagePair = { min(findImageID1, findImageID2), max(findImageID1, findImageID2) };
			if (initImagePairs.find(imagePair) != initImagePairs.end())
			{
				continue;
			}
			initImagePairs.insert(imagePair);
			if (EstimateInitialTwoViewGeometry(options, findImageID1, findImageID2))
			{
				imageID1 = findImageID1;
				imageID2 = findImageID2;
				return true;
			}
		}
	}
	return false;
}
bool CMapper::RegisterInitialImagePair(const COptions& options, size_t imageID1, size_t imageID2)
{
	Check(model.GetRegImagesNum() == 0);
	Check(database);
	options.CheckOptions();

	if (imageID1 > imageID2)
	{
		return RegisterInitialImagePair(options, imageID2, imageID1);
	}

	const pair<size_t, size_t> imagePair = { imageID1, imageID2 };
	initImagePairs.insert(imagePair);

	if (!EstimateInitialTwoViewGeometry(options, imageID1, imageID2))
	{
		return false;
	}

	CImage& image1 = database->GetImage(imageID1);
	const CCamera& camera1 = database->GetCamera(image1.GetCameraID());

	CImage& image2 = database->GetImage(imageID2);
	const CCamera& camera2 = database->GetCamera(image2.GetCameraID());

	CRigid3D& rigid3D1 = image1.GetWorldToCamera(modelID);
	CRigid3D& rigid3D2 = image2.GetWorldToCamera(modelID);
	const Eigen::Matrix3x4d worldToCamera1 = rigid3D1.ToMatrix();
	const Eigen::Matrix3x4d worldToCamera2 = rigid3D2.ToMatrix();
	const Eigen::Vector3d projectCenter1 = image1.GetProjectionCenter(modelID);
	const Eigen::Vector3d projectCenter2 = image2.GetProjectionCenter(modelID);

	rigid3D1 = CRigid3D();
	rigid3D2 = preInitTwoViewGeometry.image1ToImage2;

	model.RegisterImage(imageID1);
	model.RegisterImage(imageID2);
	RegisterImageEvent(imageID1);
	RegisterImageEvent(imageID2);

	const CSIFTMatches& matches = database->GetTwoViewGeometry(imageID1, imageID2).inlierMatches;
	const double minTriAngleRad = DegToRad(options.reconstructionOptions.initMinTriAngle);

	CTrack track;
	track.Reserve(2);
	track.AddElement(CTrackElement());
	track.AddElement(CTrackElement());
	track.GetElement(0).imageID = imageID1;
	track.GetElement(1).imageID = imageID2;

	for (const CSIFTMatch& match : matches)
	{
		const CKeypoint& keypoint1 = image1.GetKeypoint(match.point2DIndex1);
		const CKeypoint& keypoint2 = image2.GetKeypoint(match.point2DIndex2);
		const Eigen::Vector2d point1 = camera1.ImageToCamera(Eigen::Vector2d(keypoint1.pt.x, keypoint1.pt.y));
		const Eigen::Vector2d point2 = camera2.ImageToCamera(Eigen::Vector2d(keypoint2.pt.x, keypoint2.pt.y));
		const Eigen::Vector3d XYZ = TriangulatePoint(worldToCamera1, worldToCamera2, point1, point2);
		const double TriAngle = CalculateTriangulationAngle(projectCenter1, projectCenter2, XYZ);
		if (TriAngle >= minTriAngleRad && HasPointPositiveDepth(worldToCamera1, XYZ) && HasPointPositiveDepth(worldToCamera2, XYZ))
		{
			track.GetElement(0).point2DIndex = match.point2DIndex1;
			track.GetElement(1).point2DIndex = match.point2DIndex2;
			model.AddPoint3D(XYZ, track);
		}
	}
	return true;
}
vector<size_t> CMapper::FindNextImages(const COptions& options)
{
	Check(database);
	options.CheckOptions();

	const size_t numImages = database->GetImagesNum();
	vector<pair<size_t, float>> imageRanks;
	vector<pair<size_t, float>> otherImageRanks;
	imageRanks.reserve(numImages);
	otherImageRanks.reserve(numImages);
	for (size_t imageID = 0; imageID < numImages; imageID++)
	{
		if (model.IsImageRegistered(imageID))
		{
			continue;
		}
		const CImage& image = database->GetImage(imageID);
		Check(!image.IsRegistered(modelID));
		if (image.GetNumVisiblePoints3D(modelID) < options.reconstructionOptions.absPoseMinNumInliers)
		{
			continue;
		}

		float rank;
		switch (options.reconstructionOptions.nextImageSelectionMethod)
		{
		case CNextImageSelectionMethod::CMaxVisiblePointsNum:
			rank = image.GetNumVisiblePoints3D(modelID);
			break;
		case CNextImageSelectionMethod::CMaxVisiblePointsRatio:
			rank = image.GetNumVisiblePoints3D(modelID) * 1.0 / image.GetNumObservations();
			break;
		case CNextImageSelectionMethod::CMinUncertainty:
			rank = image.GetPoint3DVisibilityScore();
			break;
		}

		if (filteredImages.find(imageID) == filteredImages.end())
		{
			imageRanks.push_back({ imageID,rank });
		}
		else
		{
			otherImageRanks.push_back({ imageID,rank });
		}
	}
	vector<size_t> rankedImagesIDs;
	SortAndAppendNextImages(imageRanks, rankedImagesIDs);
	SortAndAppendNextImages(otherImageRanks, rankedImagesIDs);
	return rankedImagesIDs;
}
bool CMapper::RegisterNextImage(const COptions& options, size_t imageID)
{
	Check(model.GetRegImagesNum() >= 2);
	options.CheckOptions();
	Check(database);
	Check(!model.IsImageRegistered(imageID));

	CImage& image = database->GetImage(imageID);
	CCamera& camera = database->GetCamera(image.GetCameraID());
	Check(!image.IsRegistered(modelID));

	/*if (image.GetNumVisiblePoints3D(modelID) < options.reconstructionOptions.absPoseMinNumInliers)
	{
		return false;
	}*/

	vector<pair<size_t, size_t>> triCorrs;
	vector<Eigen::Vector2d> triPoints2D;
	vector<Eigen::Vector3d> triPoints3D;

	for (size_t point2DID = 0; point2DID < image.GetNumPoints2D(); point2DID++)
	{
		const CKeypoint& point2D = image.GetKeypoint(point2DID);

		size_t corrPoint3DID = numeric_limits<size_t>::max();
		const pair<CConjugatePoints, CObjectPoints>& correspondences = image.GetCorrespondences(point2DID);
		for (const auto& pair : correspondences.first)
		{
			if (!model.IsImageRegistered(pair.first))
			{
				continue;
			}

			const CImage& corrImage = database->GetImage(pair.first);
			if (!corrImage.IsPoint2DHasPoint3D(modelID, pair.second))
			{
				continue;
			}
			const size_t point3DID = corrImage.GetPoint3DID(pair.second, modelID);
			if (corrPoint3DID == numeric_limits<size_t>::max())
			{
				corrPoint3DID = point3DID;
			}
			else
			{
				Check(corrPoint3DID == point3DID);
			}

			const CCamera& corrCamera = database->GetCamera(corrImage.GetCameraID());
			if (corrCamera.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam))
			{
				continue;
			}

			const CPoint3D& point3D = model.GetPoint3D(point3DID);
			triCorrs.push_back({ point2DID, point3DID });
			triPoints2D.push_back(Eigen::Vector2d(point2D.pt.x, point2D.pt.y));
			triPoints3D.push_back(point3D.GetXYZ());
		}
	}

	if (triPoints2D.size() < options.reconstructionOptions.absPoseMinNumInliers)
	{
		return false;
	}

	
	COptions absPoseEstimationOptions = options;
	absPoseEstimationOptions.reconstructionOptions.absPoseRANSACOoptions.maxError = absPoseEstimationOptions.reconstructionOptions.absPoseMaxError;
	absPoseEstimationOptions.reconstructionOptions.absPoseRANSACOoptions.minInlierRatio = absPoseEstimationOptions.reconstructionOptions.absPoseMinInlierRatio;
	absPoseEstimationOptions.reconstructionOptions.absPoseRANSACOoptions.minNumTrials = 100;
	absPoseEstimationOptions.reconstructionOptions.absPoseRANSACOoptions.maxNumTrials = 10000;
	absPoseEstimationOptions.reconstructionOptions.absPoseRANSACOoptions.confidence = 0.99999; // ʹ�ø����Ŷ�������P3P RANSAC�Ĺ�����ֹ


	COptions absPoseRefinementOptions = options;
	if (numRegImagesPerCamera[image.GetCameraID()] > 0) // ����Ѿ�����һ��ʹ����ͬ�����ͼ���ע������б�����
	{
		if (camera.IsBogusParams(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam))
		{
			// ֮ǰ�Ż�����������в����߼��Ĳ���, ������ò������������¹���
			const size_t width = camera.GetWidth();
			const size_t height = camera.GetHeight();
			const double focalLength = 1.2 * max(width, height);
			camera.SetFocalLengthPrior(false);
			camera.SetParams({ focalLength, width / 2.0, height / 2.0, 0 });
			absPoseEstimationOptions.reconstructionOptions.isEstimateAbsPoseFocalLength = true;
			absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseFocalLength = true;
			absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseExtraParams = true;
		}
		else
		{
			absPoseEstimationOptions.reconstructionOptions.isEstimateAbsPoseFocalLength = false;
			absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseFocalLength = false;
			absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseExtraParams = false;
		}
	}
	else // ���֮ǰδ���Ż�. ע��, �������������֮ǰ�Ѿ������Ĺ�, ��Ӱ�񱻹��˵���, ���������ʽ����������������������¹�������
	{
		camera.SetParams(database->GetCamera(image.GetCameraID()).GetParams());
		absPoseEstimationOptions.reconstructionOptions.isEstimateAbsPoseFocalLength = !camera.IsFocalLengthPrior();
		absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseFocalLength = true;
		absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseExtraParams = true;
	}

	if (!options.reconstructionOptions.isRefineAbsPoseFocalLength)
	{
		absPoseEstimationOptions.reconstructionOptions.isEstimateAbsPoseFocalLength = false;
		absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseFocalLength = false;
	}
	if (!options.reconstructionOptions.isRefineAbsPoseExtraParams)
	{
		absPoseRefinementOptions.reconstructionOptions.isRefineAbsPoseExtraParams = false;
	}

	size_t numInliers = 0;
	vector<char> inlierMask;
	if (!EstimateAbsolutePose(absPoseEstimationOptions, triPoints2D, triPoints3D, image.GetWorldToCamera(modelID), camera, numInliers, inlierMask))
	{
		return false;
	}

	if (numInliers < options.reconstructionOptions.absPoseMinNumInliers)
	{
		return false;
	}

	if (!RefineAbsolutePose(absPoseRefinementOptions, inlierMask, triPoints2D, triPoints3D, image.GetWorldToCamera(modelID), camera))
	{
		return false;
	}

	model.RegisterImage(imageID);
	RegisterImageEvent(imageID);

	for (size_t i = 0; i < inlierMask.size(); i++)
	{
		if (inlierMask[i])
		{
			const size_t point2DID = triCorrs[i].first;
			if (!image.IsPoint2DHasPoint3D(modelID, point2DID))
			{
				const size_t point3DID = triCorrs[i].second;
				const CTrackElement trackElement(imageID, point2DID);
				model.AddObservation(point3DID, trackElement);
				triangulator->AddModifiedPoint3D(point3DID);
			}
		}
	}
	return true;
}
size_t CMapper::TriangulateImage(const COptions& options, size_t imageID)
{
	Check(triangulator);
	triangulator->TriangulateImage(options, imageID);
}
size_t CMapper::Retriangulate(const COptions& options)
{
	Check(triangulator);
	triangulator->Retriangulate(options);
}
size_t CMapper::CompleteTracks(const COptions& options)
{
	Check(triangulator);
	triangulator->CompleteAllTracks(options);
}
size_t CMapper::MergeTracks(const COptions& options)
{
	Check(triangulator);
	triangulator->MergeAllTracks(options);
}
CLocalBundleAdjustmentReport CMapper::LocalBundleAdjust(const COptions& options, size_t imageID, const unordered_set<size_t>& point3DIDs)
{
	options.CheckOptions();
	Check(database);
	Check(triangulator);

	CLocalBundleAdjustmentReport report;
	const vector<size_t> localBundle = FindLocalBundle(options, imageID);
	if (!localBundle.empty())
	{
		CBundleAdjustmentConfig bundleAdjustmentConfig;
		bundleAdjustmentConfig.AddImage(imageID);
		for (size_t localImageID : localBundle)
		{
			bundleAdjustmentConfig.AddImage(localImageID);
		}
		if (options.reconstructionOptions.isFixingExistingImages)
		{
			for (size_t localImageID : localBundle)
			{
				if (existingImages.find(localImageID) != existingImages.end())
				{
					bundleAdjustmentConfig.SetConstantCameraPose(localImageID);
				}
			}
		}

		// ȷ����Щ�����Ҫ�̶�, ���ҽ���������Ĳ���������ע���Ӱ���ڵ�ǰ�ľֲ�ƽ����
		unordered_map<size_t, size_t> numImagesPerCamera;
		for (size_t imageID : bundleAdjustmentConfig.GetAllImages())
		{
			const CImage& image = database->GetImage(imageID);
			numImagesPerCamera[image.GetCameraID()]++;
		}
		for (const auto& pair : numImagesPerCamera)
		{
			const size_t cameraID = pair.first;

			const auto it = numRegImagesPerCamera.find(cameraID);
			Check(it != numRegImagesPerCamera.end());
			const size_t numRegImages = it->second;

			if (pair.second < numRegImages)
			{
				bundleAdjustmentConfig.SetConstantCameraIntrinsics(cameraID);
			}
		}

		// �̶�7�����ɶ��Ա����ھֲ�ƽ���г��ֳ߶�/��ת/ƽ�Ƶ�Ư��
		if (localBundle.size() == 1)
		{
			bundleAdjustmentConfig.SetConstantCameraPose(localBundle[0]);
			bundleAdjustmentConfig.SetConstantCameraPositions(imageID, { 0 });
		}
		else
		{
			const size_t imageID1 = localBundle[localBundle.size() - 1];
			const size_t imageID2 = localBundle[localBundle.size() - 2];
			bundleAdjustmentConfig.SetConstantCameraPose(imageID1);
			if (!options.reconstructionOptions.isFixingExistingImages || existingImages.find(imageID2) == existingImages.end())
			{
				bundleAdjustmentConfig.SetConstantCameraPositions(imageID2, { 0 });
			}
		}

		// ȷ�������Ż������µĺͶ̹켣��3D��, ���������Ƿ���ȫ�����ڱ���ƽ����.
		// ���������켣��3D��, ��Ϊ����ͨ�����Ѿ��ǳ��ȶ�
		unordered_set<size_t> variablePoint3DIDs;
		for (size_t point3DID : point3DIDs)
		{
			const CPoint3D& point3D = model.GetPoint3D(point3DID);
			if (!point3D.HasError() || point3D.GetTrack().GetTrackLength() <= 15)
			{
				bundleAdjustmentConfig.AddVariablePoint(point3DID);
				variablePoint3DIDs.insert(point3DID);
			}
		}

		CBundleAdjuster bundleAdjuster(options, bundleAdjustmentConfig, database);
		bundleAdjuster.Solve(model);

		report.numAdjustedObservations = bundleAdjuster.GetSummary().num_residuals / 2;
		report.numMergedObservations = triangulator->MergeTracks(options, variablePoint3DIDs);
		report.numCompletedObservations = triangulator->CompleteTracks(options, variablePoint3DIDs);
		report.numCompletedObservations += triangulator->CompleteImage(options, imageID);
	}

	// ���˱��޸ĵ�ͼ������иı��3D��, ��ȷ��ģ����û���쳣��.
	// ��ᵼ���ظ�����, ��Ϊ����ṩ��3D�����Ҳ�����ڵ�������ͼ����, ��������׶�, ���˲�����һ��ƿ��
	unordered_set<size_t> filterImageIDs;
	filterImageIDs.insert(imageID);
	filterImageIDs.insert(localBundle.begin(), localBundle.end());

	report.numFilteredObservations = model.FilterPoints3DInImages(options.reconstructionOptions.filterMaxReprojectionError, options.reconstructionOptions.filterMinTriAngle, filterImageIDs);
	report.numFilteredObservations += model.FilterPoints3D(options.reconstructionOptions.filterMaxReprojectionError, options.reconstructionOptions.filterMinTriAngle, point3DIDs);

	return report;
}
bool CMapper::GlobalBundleAdjust(const COptions& options)
{
	options.CheckOptions();
	Check(database);
	Check(triangulator);

	// ���ﲻ��ʹ��model.GetAllRegImages(), ��Ϊ����Ҫ���մ���Ӱ���˳������, ������Ҫ�ѵ�0��Ӱ������λ�˹̶�
	vector<size_t> regImageIDs;
	regImageIDs.reserve(model.GetRegImagesNum());
	for (size_t imageID = 0; imageID < database->GetImagesNum(); imageID++)
	{
		if (model.IsImageRegistered(imageID))
		{
			regImageIDs.push_back(imageID);
		}
	}
	Check(regImageIDs.size() >= 2);

	model.FilterObservationsWithNegativeDepth();

	CBundleAdjustmentConfig bundleAdjustmentConfig;
	for (size_t imageID : regImageIDs)
	{
		bundleAdjustmentConfig.AddImage(imageID);
	}

	if (options.reconstructionOptions.isFixingExistingImages)
	{
		for (size_t imageID : regImageIDs)
		{
			if (existingImages.find(imageID) != existingImages.end())
			{
				bundleAdjustmentConfig.SetConstantCameraPose(imageID);
			}
		}
	}

	bundleAdjustmentConfig.SetConstantCameraPose(regImageIDs[0]);
	if (!options.reconstructionOptions.isFixingExistingImages || existingImages.find(regImageIDs[1]) == existingImages.end())
	{
		bundleAdjustmentConfig.SetConstantCameraPositions(regImageIDs[1], { 0 });
	}

	CBundleAdjuster bundleAdjuster(options, bundleAdjustmentConfig, database);
	if (!bundleAdjuster.Solve(model))
	{
		return false;
	}

	// �Գ������й�һ���������ֵ�ȶ���, ���ұ�����ģ��������г��ִ�ĳ߶ȱ仯
	model.Normalize();

	return true;
}
size_t CMapper::FilterImages(const COptions& options)
{
	options.CheckOptions();

	// ���ؽ������ڽ׶β�Ҫ����ͼ��, ��Ϊ������ڲ�ͨ�����ڴ���Ż���. ���, ��������ڿ�ʼʱ�����ȶ�
	if (model.GetRegImagesNum() < 20)
	{
		return 0;
	}

	const vector<size_t> imageIDs = model.FilterImages(options.reconstructionOptions.minFocalLengthRatio, options.reconstructionOptions.maxFocalLengthRatio, options.reconstructionOptions.maxExtraParam);
	for (size_t imageID : imageIDs)
	{
		DeRegisterImageEvent(imageID);
		filteredImages.insert(imageID);
	}
	return imageIDs.size();
}
size_t CMapper::FilterPoints(const COptions& options)
{
	options.CheckOptions();
	return model.FilterAllPoints3D(options.reconstructionOptions.filterMaxReprojectionError, options.reconstructionOptions.filterMinTriAngle);
}
const CModel& CMapper::GetModel() const
{
	return model;
}
const unordered_set<size_t>& CMapper::GetModifiedPoints3D()
{
	Check(triangulator);
	return triangulator->GetModifiedPoints3D();
}
void CMapper::ClearModifiedPoints3D()
{
	Check(triangulator);
	triangulator->ClearModifiedPoints3D();
}
vector<size_t> CMapper::FindFirstInitialImage() const
{
	Check(database);
	struct CImageInfo
	{
		size_t imageID;
		bool isFocalLengthPrior;
		size_t numCorrespondences;
	};

	vector<CImageInfo> imageInfos;
	imageInfos.reserve(database->GetImagesNum());
	for (size_t imageID = 0; imageID < database->GetImagesNum(); imageID++)
	{
		const CImage& image = database->GetImage(imageID);
		if (image.GetNumCorrespondences() == 0)
		{
			continue;
		}
		if (image.GetNumRegisteredModels() > 0)
		{
			continue;
		}
		const CCamera& camera = database->GetCamera(image.GetCameraID());

		CImageInfo imageInfo;
		imageInfo.imageID = imageID;
		imageInfo.isFocalLengthPrior = camera.IsFocalLengthPrior();
		imageInfo.numCorrespondences = image.GetNumCorrespondences();
		imageInfos.push_back(imageInfo);
	}

	sort(imageInfos.begin(), imageInfos.end(), [](const CImageInfo& imageInfo1, const CImageInfo& imageInfo2)
		{
			if (imageInfo1.isFocalLengthPrior && !imageInfo2.isFocalLengthPrior)
			{
				return true;
			}
			if (!imageInfo1.isFocalLengthPrior && imageInfo2.isFocalLengthPrior)
			{
				return true;
			}
			return imageInfo1.numCorrespondences > imageInfo2.numCorrespondences;
		});

	vector<size_t> imageIDs(imageInfos.size());
	for (size_t i = 0; i < imageInfos.size(); i++)
	{
		imageIDs[i] = imageInfos[i].imageID;
	}
	return imageIDs;
}
vector<size_t> CMapper::FindSecondInitialImage(const COptions& options, size_t imageID1) const
{
	options.CheckOptions();
	Check(database);
	struct CImageInfo
	{
		size_t imageID;
		bool isFocalLengthPrior;
		size_t numCorrespondences;
	};

	const CImage& image1 = database->GetImage(imageID1);
	const size_t numImages = database->GetImagesNum();

	// �ռ����һ�ų�ʼӰ�����Ӳ���û�б��κ�ģ��ע�����ͼ��
	vector<CImageInfo> imageInfos;
	imageInfos.reserve(numImages);
	for (size_t imageID = 0; imageID < numImages; imageID++)
	{
		if (imageID == imageID1 || database->GetImage(imageID).GetNumRegisteredModels() > 0)
		{
			continue;
		}
		const size_t inlierMatchesNum = database->GetInlierMatchesNum(imageID, imageID1);
		if (inlierMatchesNum >= options.reconstructionOptions.initMinNumInliers)
		{
			const CImage& image = database->GetImage(imageID);
			const CCamera& camera = database->GetCamera(image.GetCameraID());

			CImageInfo imageInfo;
			imageInfo.imageID = imageID;
			imageInfo.isFocalLengthPrior = camera.IsFocalLengthPrior();
			imageInfo.numCorrespondences = inlierMatchesNum;
			imageInfos.push_back(imageInfo);
		}
	}

	sort(imageInfos.begin(), imageInfos.end(), [](const CImageInfo& imageInfo1, const CImageInfo& imageInfo2)
		{
			if (imageInfo1.isFocalLengthPrior && !imageInfo2.isFocalLengthPrior)
			{
				return true;
			}
			if (!imageInfo1.isFocalLengthPrior && imageInfo2.isFocalLengthPrior)
			{
				return true;
			}
			return imageInfo1.numCorrespondences > imageInfo2.numCorrespondences;
		});

	vector<size_t> imageIDs(imageInfos.size());
	for (size_t i = 0; i < imageInfos.size(); i++)
	{
		imageIDs[i] = imageInfos[i].imageID;
	}
	return imageIDs;
}
vector<size_t> CMapper::FindLocalBundle(const COptions& options, size_t imageID) const
{
	options.CheckOptions();
	Check(model.IsImageRegistered(imageID));
	Check(database);

	const CImage& image = database->GetImage(imageID);

	// �������Ӱ��������һ����ͬ3D���Ӱ��, ���Ҽ��㹲ͬ3D�������
	unordered_map<size_t, size_t> sharedObservations;
	sharedObservations.reserve(database->GetImagesNum());
	unordered_set<size_t> point3DIDs;
	point3DIDs.reserve(model.GetPoints3DNum());
	for (size_t point2DID = 0; point2DID < image.GetNumPoints2D(); point2DID++)
	{
		if (!image.IsPoint2DHasPoint3D(modelID, point2DID))
		{
			continue;
		}
		const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
		const CPoint3D& point3D = model.GetPoint3D(point3DID);
		const vector<CTrackElement>& trackElements = point3D.GetTrack().GetAllElements();
		for (const CTrackElement& trackElement : trackElements)
		{
			if (trackElement.imageID != imageID)
			{
				sharedObservations[trackElement.imageID]++;
			}
		}
	}

	// ���ݹ���۲����������ص�ͼ���������
	vector<pair<size_t, size_t>> overlappingImages(sharedObservations.begin(), sharedObservations.end());
	sort(overlappingImages.begin(), overlappingImages.end(), [](const pair<size_t, size_t>& image1, const pair<size_t, size_t>& image2)
		{
			return image1.second > image2.second;
		});

	// ����Ӱ���ɸ���ͼ��������ص��ھ�ͼ�����, ��˼�ȥ1
	const size_t numImages = options.reconstructionOptions.numLocalBundleAdjustmentImages - 1;
	const size_t numEffImages = min(numImages, overlappingImages.size());

	// ��ȡ���������ͼ��ȷ�����㹻�����ǻ��Ƕ�
	vector<size_t> localBundleImageIDs;
	localBundleImageIDs.reserve(numEffImages);
	if (overlappingImages.size() == numEffImages)
	{
		for (const pair<size_t, size_t>& pair : overlappingImages)
		{
			localBundleImageIDs.push_back(pair.first);
		}
		return localBundleImageIDs;
	}

	// �ڽ������ĵ�����, ���Ǵ��ص�����ߵ�ͼ��ʼ, ������Ƿ�����㹻�����ǻ��Ƕ�. ���û��һ���ص�ͼ������㹻�����ǻ��Ƕ�, ���ǽ��ſ����ǻ��Ƕȵ���ֵ, ���ٴδ��ص�����ߵ�ͼ��ʼ
	// ���գ����������Ȼû���ҵ��㹻��ͼ��, ���Ǿͼ򵥵�ʹ���ص�����ߵ�ͼ��
	const double minTriAngleRad = DegToRad(options.reconstructionOptions.minTriAngleLocalBundleAdjustment);

	// ���ηſ��ѡ����ֵ(��С���ǻ��Ƕ�, ��С����۲������)
	const array<pair<double, double>, 8> selectionThresholds = { {
		{minTriAngleRad / 1.0, 0.6 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 1.5, 0.6 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 2.0, 0.5 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 2.5, 0.4 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 3.0, 0.3 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 4.0, 0.2 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 5.0, 0.1 * image.GetNumPoints3D(modelID)},
		{minTriAngleRad / 6.0, 0.1 * image.GetNumPoints3D(modelID)},
		} };

	const Eigen::Vector3d projectionCenter = image.GetProjectionCenter(modelID);
	vector<Eigen::Vector3d> sharedPoints3D;
	sharedPoints3D.reserve(image.GetNumPoints3D(modelID));
	vector<double> triAngles(overlappingImages.size(), -1);
	vector<char> usedOverlappingImages(overlappingImages.size(), false);

	for (const pair<double, double>& selectedThreshold : selectionThresholds)
	{
		for (size_t overlappingImageIndex = 0; overlappingImageIndex < overlappingImages.size(); overlappingImageIndex++)
		{
			const size_t overlappingImageID = overlappingImages[overlappingImageIndex].first;
			const size_t numCommonPoints3D = overlappingImages[overlappingImageIndex].second;
			// ���ͼ���Ƿ����㹻���ص�. ����ͼ���ǻ����ص����������, ���ǿ���ֱ������ʣ���ͼ��
			if (numCommonPoints3D < selectedThreshold.second)
			{
				break;
			}
			if (usedOverlappingImages[overlappingImageIndex])
			{
				continue;
			}
			const CImage& overlappingImage = database->GetImage(overlappingImageID);
			const Eigen::Vector3d overlappingImageProjectionCenter = overlappingImage.GetProjectionCenter(modelID);

			// �ڵ�һ�ε�����, �������ǻ��Ƕ�. �ں����ĵ�����, ������ǰ�����ֵ
			if (triAngles[overlappingImageIndex] < 0)
			{
				// �ռ��ܹ�ͬ�۲⵽��3D��
				sharedPoints3D.clear();
				const size_t numOverlappingImagePoints2D = overlappingImage.GetNumPoints2D();
				for (size_t point2DID = 0; point2DID < numOverlappingImagePoints2D; point2DID++)
				{
					if (!overlappingImage.IsPoint2DHasPoint3D(modelID, point2DID))
					{
						continue;
					}
					const size_t point3DID = overlappingImage.GetPoint3DID(point2DID, modelID);
					if (point3DIDs.find(point3DID) != point3DIDs.end())
					{
						sharedPoints3D.push_back(model.GetPoint3D(point3DID).GetXYZ());
					}
				}
				// ����75%��λ�����Ľ����
				triAngles[overlappingImageIndex] = Percentile(CalculateTriangulationAngles(projectionCenter, overlappingImageProjectionCenter, sharedPoints3D), 75);
			}
			
			// ���Ӱ���Ƿ�����㹻�Ľ����
			if (triAngles[overlappingImageIndex] >= selectedThreshold.first)
			{
				localBundleImageIDs.push_back(overlappingImageID);
				usedOverlappingImages[overlappingImageIndex] = true;

				// ����Ƿ��Ѿ��ռ����㹻��ͼ��
				if (localBundleImageIDs.size() >= numEffImages)
				{
					break;
				}
			}
		}
		if (localBundleImageIDs.size() >= numEffImages)
		{
			break;
		}
	}

	// ���û���㹻��ľ����㹻��Ľ���ǵ�Ӱ��, ��򵥵�������ص���Ӱ�����ʣ�ಿ��
	if (localBundleImageIDs.size() < numEffImages)
	{
		for (size_t overlappingImageIndex = 0; overlappingImageIndex < overlappingImages.size(); overlappingImageIndex++)
		{
			if (usedOverlappingImages[overlappingImageIndex])
			{
				continue;
			}

			// ���Ӱ����δ�ھֲ�ƽ�Χ��, ���ռ���ͼ��
			localBundleImageIDs.push_back(overlappingImages[overlappingImageIndex].first);
			usedOverlappingImages[overlappingImageIndex] = true;

			if (localBundleImageIDs.size() >= numEffImages)
			{
				break;
			}
		}
	}
	return localBundleImageIDs;
}
void CMapper::RegisterImageEvent(size_t imageID)
{
	Check(database);
	const CImage& image = database->GetImage(imageID);
	numRegImagesPerCamera[image.GetCameraID()]++;
}
void CMapper::DeRegisterImageEvent(size_t imageID)
{
	Check(database);
	const CImage& image = database->GetImage(imageID);

	auto it = numRegImagesPerCamera.find(image.GetCameraID());
	Check(it != numRegImagesPerCamera.end() && it->second > 0);
	it->second--;
}
bool CMapper::EstimateInitialTwoViewGeometry(const COptions& options, size_t imageID1, size_t imageID2)
{
	if (imageID1 > imageID2)
	{
		return EstimateInitialTwoViewGeometry(options, imageID2, imageID1);
	}
	if (preInitImagePair.first == imageID1 && preInitImagePair.second == imageID2)
	{
		return true;
	}
	const CImage& image1 = database->GetImage(imageID1);
	const CImage& image2 = database->GetImage(imageID2);

	const CCamera& camera1 = database->GetCamera(image1.GetCameraID());
	const CCamera& camera2 = database->GetCamera(image2.GetCameraID());

	














}





















