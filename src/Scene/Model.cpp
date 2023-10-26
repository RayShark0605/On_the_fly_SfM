#include "Model.h"

using namespace std;


CModel::CModel(size_t modelID, CDatabase* database, const COptions& options) :database(database)
{
	Check(database);
	options.CheckOptions();

	this->modelID = modelID;
	regImageIDs.clear();
	points3D.clear();
	nextPoint3DID = 0;

	const unordered_map<pair<size_t, size_t>, size_t, MatchPairHash, MatchPairEqual> correspondences = database->GetAllCorrespondences();
	for (const auto& pair : correspondences)
	{
		const size_t imageID1 = pair.first.first;
		const size_t imageID2 = pair.first.second;
		const size_t numCorrespondences = pair.second;
		if (numCorrespondences >= options.reconstructionOptions.minNumMatches)
		{
			Check(imageID1 < imageID2);
			CImagePairStatus imagePairStatus;
			imagePairStatus.numTotalCorrs = numCorrespondences;
			imagePairs[pair.first] = imagePairStatus;
		}
	}
}
unordered_set<size_t> CModel::GetAllPoints3D() const
{
	unordered_set<size_t> result;
	result.reserve(points3D.size());
	for (const auto& pair : points3D)
	{
		result.insert(pair.first);
	}
	return result;
}
size_t CModel::AddPoint3D(const Eigen::Vector3d& XYZ, const CTrack& track, const Eigen::Vector3ub& color)
{
	Check(database);
	size_t point3DID = nextPoint3DID;
	while (points3D.find(nextPoint3DID) != points3D.end())
	{
		nextPoint3DID++;
	}
	Check(point3DID != nextPoint3DID);

	const vector<CTrackElement> trackElements = track.GetAllElements();
	for (const CTrackElement& trackElement : trackElements)
	{
		Check(regImageIDs.find(trackElement.imageID) != regImageIDs.end());
		CImage& image = database->GetImage(trackElement.imageID);
		Check(!image.IsPoint2DHasPoint3D(modelID, trackElement.point2DIndex));
		image.SetPoint3DForPoint2D(trackElement.point2DIndex, point3DID, modelID);
		Check(image.GetNumPoints2D() >= image.GetNumPoints3D(modelID));
	}
	for (const CTrackElement& trackElement : trackElements)
	{
		SetObservationAsTriangulated(trackElement.imageID, trackElement.point2DIndex, false);
	}

	points3D[point3DID] = CPoint3D(XYZ);
	points3D[point3DID].SetTrack(track);
	points3D[point3DID].SetColor(color);
	return point3DID;
}
void CModel::AddObservation(size_t point3DID, const CTrackElement& trackElement)
{
	Check(database);
	Check(points3D.find(point3DID) != points3D.end());

	CImage& image = database->GetImage(trackElement.imageID);
	Check(!image.IsPoint2DHasPoint3D(modelID, trackElement.point2DIndex));
	image.SetPoint3DForPoint2D(trackElement.point2DIndex, point3DID, modelID);
	Check(image.GetNumPoints3D(modelID) <= image.GetNumPoints2D());

	points3D[point3DID].GetTrack().AddElement(trackElement);
	SetObservationAsTriangulated(trackElement.imageID, trackElement.point2DIndex, true);
}
size_t CModel::MergePoints3D(size_t point3D1ID, size_t point3D2ID)
{
	Check(points3D.find(point3D1ID) != points3D.end() && points3D.find(point3D2ID) != points3D.end());

	const CPoint3D& point1 = points3D[point3D1ID];
	const CPoint3D& point2 = points3D[point3D2ID];

	const size_t trackLength1 = point1.GetTrack().GetTrackLength();
	const size_t trackLength2 = point2.GetTrack().GetTrackLength();

	const Eigen::Vector3d mergedXYZ = (trackLength1 * point1.GetXYZ() + trackLength2 * point2.GetXYZ()) / (trackLength1 + trackLength2);
	const Eigen::Vector3d mergedRGB = (trackLength1 * point1.GetColor().cast<double>() + trackLength2 * point2.GetColor().cast<double>()) / (trackLength1 + trackLength2);

	CTrack mergedTrack;
	mergedTrack.Reserve(trackLength1 + trackLength2);
	mergedTrack.AddElements(point1.GetTrack().GetAllElements());
	mergedTrack.AddElements(point2.GetTrack().GetAllElements());
	DeletePoint3D(point3D1ID);
	DeletePoint3D(point3D2ID);
	const size_t mergedPoint3DID = AddPoint3D(mergedXYZ, mergedTrack, mergedRGB.cast<uint8_t>());
	return mergedPoint3DID;
}
void CModel::DeletePoint3D(size_t point3DID)
{
	// 注意：不要更改这些代码的顺序
	Check(database);
	Check(points3D.find(point3DID) != points3D.end());
	const CTrack& track = points3D[point3DID].GetTrack();
	const vector<CTrackElement> trackElements = track.GetAllElements();
	for (const CTrackElement& trackElement : trackElements)
	{
		ResetTriObservations(trackElement.imageID, trackElement.point2DIndex, true);
	}
	for (const CTrackElement& trackElement : trackElements)
	{
		database->GetImage(trackElement.imageID).ResetPoint3DForPoint2D(trackElement.point2DIndex, modelID);
	}
	points3D.erase(point3DID);
	nextPoint3DID = min(nextPoint3DID, point3DID);
}
void CModel::DeleteObservation(size_t imageID, size_t point2DID)
{
	// 注意：不要更改这些代码的顺序
	Check(database);
	CImage& image = database->GetImage(imageID);
	const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
	
	const auto it = points3D.find(point3DID);
	Check(it != points3D.end());
	CPoint3D& point3D = it->second;
	CTrack& track = point3D.GetTrack();
	if (track.GetTrackLength() <= 2)
	{
		DeletePoint3D(point3DID);
		return;
	}

	track.DeleteElement(imageID, point2DID);
	ResetTriObservations(imageID, point2DID, false);
	image.ResetPoint3DForPoint2D(point2DID, modelID);
}
void CModel::RegisterImage(size_t imageID)
{
	Check(database);
	CImage& image = database->GetImage(imageID);
	const auto it = regImageIDs.find(modelID);
	if (!image.IsRegistered(modelID))
	{
		Check(it == regImageIDs.end());
		image.SetRegistered(modelID, true);
		regImageIDs.insert(imageID);
	}
	Check(it != regImageIDs.end());
}
void CModel::DeRegisterImage(size_t imageID)
{
	Check(database);
	CImage& image = database->GetImage(imageID);
	const size_t numPoints2D = image.GetNumPoints2D();
	for (size_t point2DID = 0; point2DID < numPoints2D; point2DID++)
	{
		if (image.IsPoint2DHasPoint3D(modelID, point2DID))
		{
			DeleteObservation(imageID, point2DID);
		}
	}
	image.SetRegistered(modelID, false);
	regImageIDs.erase(imageID);
}
void CModel::Normalize(double extent, double p0, double p1, bool isUseImage)
{
	if ((isUseImage && regImageIDs.size() < 2) || (!isUseImage && points3D.size() < 2))
	{
		return;
	}
	Check(extent > 0 && p0 >= 0 && p0 <= 1 && p1 >= 0 && p1 <= 1 && p0 <= p1);
	const tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> boundsAndCentroid = ComputeBoundsAndCentroid(p0, p1, isUseImage);

	// 计算平移和缩放量, 以便在缩放之前执行平移
	const double oldExtent = (get<1>(boundsAndCentroid) - get<0>(boundsAndCentroid)).norm();
	double scale;
	if (oldExtent < numeric_limits<double>::epsilon())
	{
		scale = 1;
	}
	else
	{
		scale = extent / oldExtent;
	}
	CSim3D tform(scale, Eigen::Quaterniond::Identity(), -scale * get<2>(boundsAndCentroid));
	Transform(tform);
}
Eigen::Vector3d CModel::ComputeCentroid(double p0, double p1) const
{
	Check(p0 >= 0 && p0 <= 1 && p1 >= 0 && p1 <= 1 && p0 <= p1);
	return get<2>(ComputeBoundsAndCentroid(p0, p1, false));
}
pair<Eigen::Vector3d, Eigen::Vector3d> CModel::ComputeBoundingBox(double p0, double p1) const
{
	Check(p0 >= 0 && p0 <= 1 && p1 >= 0 && p1 <= 1 && p0 <= p1);
	const tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> boundsAndCentroid = ComputeBoundsAndCentroid(p0, p1, false);
	return { get<0>(boundsAndCentroid), get<1>(boundsAndCentroid) };
}
void CModel::Transform(const CSim3D& oldWorldToNewWorld)
{
	Check(database);
	// TODO: 原Colmap的Transform函数会转换整个database所有的image, 本人认为仅需要转换当前模型注册的影像即可, 不被当前影像注册的影像应该不用转换
	for (size_t regImageID : regImageIDs)
	{
		CImage& image = database->GetImage(regImageID);
		image.GetWorldToCamera(modelID) = TransformCameraWorld(oldWorldToNewWorld, image.GetWorldToCamera(modelID));
	}
	for (auto& pair : points3D)
	{
		pair.second.SetXYZ(oldWorldToNewWorld * pair.second.GetXYZ());
	}
}
vector<size_t> CModel::FindCommonRegisteredImages(const CModel& otherModel) const
{
	Check(database && database == otherModel.database);
	vector<size_t> result;
	result.reserve(regImageIDs.size());
	for (size_t regImageID : regImageIDs)
	{
		const CImage& image = database->GetImage(regImageID);
		if (otherModel.regImageIDs.find(regImageID) != otherModel.regImageIDs.end())
		{
			Check(image.IsRegistered(modelID) && image.IsRegistered(otherModel.modelID));
			result.push_back(regImageID);
		}
		else
		{
			Check(image.IsRegistered(modelID) && !image.IsRegistered(otherModel.modelID));
		}
	}
	return result;
}
size_t CModel::FilterPoints3D(double maxReprojectionError, double minTriAngle, const unordered_set<size_t>& pointsToBeFiltered)
{
	size_t numFiltered = 0;
	numFiltered += FilterPoints3DWithLargeReprojectionError(maxReprojectionError, pointsToBeFiltered);
	numFiltered += FilterPoints3DWithSmallTriangulationAngle(minTriAngle, pointsToBeFiltered);
	return numFiltered;
}
size_t CModel::FilterPoints3DInImages(double maxReprojectionError, double minTriAngle, const unordered_set<size_t>& imageIDs)
{
	Check(database);
	unordered_set<size_t> pointsToBeFiltered;
	for (size_t imageID : imageIDs)
	{
		const CImage& image = database->GetImage(imageID);
		const size_t numPoints2D = image.GetNumPoints2D();
		for (size_t point2DID = 0; point2DID < numPoints2D; point2DID++)
		{
			if (image.IsPoint2DHasPoint3D(modelID, point2DID))
			{
				const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
				Check(points3D.find(point3DID) != points3D.end());
				pointsToBeFiltered.insert(point3DID);
			}
		}
	}
	return FilterPoints3D(maxReprojectionError, minTriAngle, pointsToBeFiltered);
}
size_t CModel::FilterAllPoints3D(double maxReprojectionError, double minTriAngle)
{
	// 首先应该过滤掉重投影误差太大的观测和点
	unordered_set<size_t> allPoints3D = GetAllPoints3D();
	size_t numFiltered = 0;
	numFiltered += FilterPoints3DWithLargeReprojectionError(maxReprojectionError, allPoints3D);
	numFiltered += FilterPoints3DWithSmallTriangulationAngle(minTriAngle, allPoints3D);
	return numFiltered;
}
size_t CModel::FilterObservationsWithNegativeDepth()
{
	Check(database);
	size_t numFiltered = 0;
	for (size_t imageID : regImageIDs)
	{
		const CImage& image = database->GetImage(imageID);
		const size_t numPoints2D = image.GetNumPoints2D();
		const Eigen::Matrix3x4d worldToCamera = image.GetWorldToCamera(modelID).ToMatrix();
		for (size_t point2DID = 0; point2DID < numPoints2D; point2DID++)
		{
			if (image.IsPoint2DHasPoint3D(modelID, point2DID))
			{
				const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
				Check(points3D.find(point3DID) != points3D.end());
				const CPoint3D& point3D = points3D[point3DID];
				if (!HasPointPositiveDepth(worldToCamera, point3D.GetXYZ()))
				{
					DeleteObservation(imageID, point2DID);
					numFiltered++;
				}
			}
		}
	}
	return numFiltered;
}
vector<size_t> CModel::FilterImages(double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam)
{
	Check(database);
	vector<size_t> filteredImageIDs;
	for (size_t imageID : regImageIDs)
	{
		const CImage& image = database->GetImage(imageID);
		if (image.GetNumPoints3D(modelID) == 0 || database->GetCamera(image.GetCameraID()).IsBogusParams(minFocalLengthRatio, maxFocalLengthRatio, maxExtraParam))
		{
			filteredImageIDs.push_back(imageID);
		}
	}
	for (size_t imageID : filteredImageIDs)
	{
		DeRegisterImage(imageID);
	}
	return filteredImageIDs;
}
size_t CModel::ComputeNumObservations() const
{
	Check(database);
	size_t numObservations = 0;
	for (size_t imageID : regImageIDs)
	{
		const CImage& image = database->GetImage(imageID);
		numObservations += image.GetNumPoints3D(modelID);
	}
	return numObservations;
}
double CModel::ComputeMeanTrackLength() const
{
	if (points3D.empty()) return 0;
	return ComputeNumObservations() * 1.0 / points3D.size();
}
double CModel::ComputeMeanObservationsPerRegImage() const
{
	if (regImageIDs.empty()) return 0;
	return ComputeNumObservations() * 1.0 / regImageIDs.size();
}
double CModel::ComputeMeanReprojectionError() const
{
	double errorSum = 0;
	size_t numValidError = 0;
	for (const auto& pair : points3D)
	{
		if (pair.second.HasError())
		{
			double error = pair.second.GetError();
			Check(error >= 0);
			errorSum += error;
			numValidError++;
		}
	}
	if (numValidError == 0)
	{
		return 0;
	}
	return errorSum / numValidError;
}
void CModel::UpdatePoint3DErrors()
{
	Check(database);
	for (auto& pair : points3D)
	{
		if (pair.second.GetTrack().GetTrackLength() == 0)
		{
			pair.second.SetError(0);
			continue;
		}
		double errorSum = 0;
		const vector<CTrackElement>& trackElements = pair.second.GetTrack().GetAllElements();
		for (const CTrackElement& trackElement : trackElements)
		{
			const CImage& image = database->GetImage(trackElement.imageID);
			const CKeypoint& point2D = image.GetKeypoint(trackElement.point2DIndex);
			const CCamera& camera = database->GetCamera(image.GetCameraID());
			const Eigen::Vector2d XY(point2D.pt.x, point2D.pt.y);
			errorSum += sqrt(CalculateSquaredReprojectionError(XY, pair.second.GetXYZ(), image.GetWorldToCamera(modelID), camera));
		}
		pair.second.SetError(errorSum / trackElements.size());
	}
}
bool CModel::ExtractColorsForImage(size_t imageID, const string& imagePath)
{
	Check(database);
	if (!IsFileExists(imagePath))
	{
		return false;
	}
	const cv::Mat imageMat = cv::imread(imagePath, cv::IMREAD_COLOR);
	if (imageMat.empty())
	{
		return false;
	}
	const Eigen::Vector3ub blackColor(0, 0, 0);

	auto bilinearInterpolation = [](const cv::Mat& img, float x, float y) -> Eigen::Vector3ub
		{
		// 获取四个临近像素的整数坐标
		int x1 = floor(x);
		int x2 = ceil(x);
		int y1 = floor(y);
		int y2 = ceil(y);

		// 确保坐标在图像范围内
		x1 = max(0, min(x1, img.cols - 1));
		x2 = max(0, min(x2, img.cols - 1));
		y1 = max(0, min(y1, img.rows - 1));
		y2 = max(0, min(y2, img.rows - 1));

		// 获取四个临近像素的值
		cv::Vec3b p1 = img.at<cv::Vec3b>(y1, x1);
		cv::Vec3b p2 = img.at<cv::Vec3b>(y1, x2);
		cv::Vec3b p3 = img.at<cv::Vec3b>(y2, x1);
		cv::Vec3b p4 = img.at<cv::Vec3b>(y2, x2);

		// 进行双线性插值
		float dx1 = x - x1;
		float dx2 = x2 - x;
		float dy1 = y - y1;
		float dy2 = y2 - y;

		cv::Vec3b interpolated = (dy2 * (dx2 * p1 + dx1 * p2) + dy1 * (dx2 * p3 + dx1 * p4));

		return Eigen::Vector3ub(interpolated[2], interpolated[1], interpolated[0]);
		};


	const CImage& image = database->GetImage(imageID);
	const size_t numPoints2D = image.GetNumPoints2D();
	for (size_t point2DID = 0; point2DID < numPoints2D; point2DID++)
	{
		if (image.IsPoint2DHasPoint3D(modelID, point2DID))
		{
			const CKeypoint& point2D = image.GetKeypoint(point2DID);
			const size_t point3DID = image.GetPoint3DID(point2DID, modelID);
			Check(points3D.find(point3DID) != points3D.end());
			CPoint3D& point3D = points3D[point3DID];
			if (point3D.GetColor() != blackColor)
			{
				continue;
			}
			point3D.SetColor(bilinearInterpolation(imageMat, point2D.pt.x, point2D.pt.y));
		}
	}
	return true;
}

void CModel::SetObservationAsTriangulated(size_t imageID, size_t point2DID, bool isContinuedPoint3D)
{
	Check(database);
	const CImage& image = database->GetImage(imageID);
	Check(image.IsRegistered(modelID));

	Check(image.IsPoint2DHasPoint3D(modelID, point2DID));
	const CKeypoint& point2D = image.GetKeypoint(point2DID);
	const pair<CConjugatePoints, CObjectPoints>& correspondences = image.GetCorrespondences(point2DID);
	const unordered_map<size_t, size_t>& conjugatePoints = correspondences.first;
	for (const auto& pair : conjugatePoints)
	{
		CImage& matchedImage = database->GetImage(pair.first);
		const CKeypoint& conjugatePoint = matchedImage.GetKeypoint(pair.second);

		matchedImage.IncrementCorrespondenceHasPoint3D(modelID, pair.second);
		
		const size_t thisPoint3DID = image.GetPoint3DID(point2DID, modelID);
		const size_t correspondencePoint3DID = matchedImage.GetPoint3DID(pair.second, modelID);
		Check(points3D.find(thisPoint3DID) != points3D.end());
		Check(points3D.find(correspondencePoint3DID) != points3D.end());

		// 更新像对之间共享的3D点的数量, 并确保只计算一次对应关系(只在imageID<pair.first时计算一次, 当imageID>pair.first时不计算)
		if (thisPoint3DID == correspondencePoint3DID && (isContinuedPoint3D || imageID < pair.first))
		{
			CImagePairStatus& imagePairStatus = imagePairs[make_pair(imageID, pair.first)];
			imagePairStatus.numTriCorrs++;
			Check(imagePairStatus.numTriCorrs <= imagePairStatus.numTotalCorrs);
		}
	}
}
void CModel::ResetTriObservations(size_t imageID, size_t point2DID, bool isDeletedPoint3D)
{
	Check(database);
	const CImage& image = database->GetImage(imageID);
	Check(image.IsRegistered(modelID));

	Check(image.IsPoint2DHasPoint3D(modelID, point2DID));
	const CKeypoint& point2D = image.GetKeypoint(point2DID);

	const pair<CConjugatePoints, CObjectPoints>& correspondences = image.GetCorrespondences(point2DID);
	const unordered_map<size_t, size_t>& conjugatePoints = correspondences.first;
	for (const auto& pair : conjugatePoints)
	{
		CImage& matchedImage = database->GetImage(pair.first);
		const CKeypoint& conjugatePoint = matchedImage.GetKeypoint(pair.second);

		matchedImage.DecrementCorrespondenceHasPoint3D(modelID, pair.second);

		const size_t thisPoint3DID = image.GetPoint3DID(point2DID, modelID);
		const size_t correspondencePoint3DID = matchedImage.GetPoint3DID(pair.second, modelID);
		Check(points3D.find(thisPoint3DID) != points3D.end());
		Check(points3D.find(correspondencePoint3DID) != points3D.end());
		// 更新像对之间共享的3D点的数量, 并确保只计算一次对应关系(只在imageID<pair.first时计算一次, 当imageID>pair.first时不计算)
		if (thisPoint3DID == correspondencePoint3DID && (!isDeletedPoint3D || imageID < pair.first))
		{
			CImagePairStatus& imagePairStatus = imagePairs[make_pair(imageID, pair.first)];
			Check(imagePairStatus.numTriCorrs > 0);
			imagePairStatus.numTriCorrs--;
		}
	}
}
size_t CModel::FilterPoints3DWithSmallTriangulationAngle(double minTriAngle, const unordered_set<size_t>& points3DID)
{
	Check(database);
	size_t numFiltered = 0;
	const double minTriAngle_Rad = DegToRad(minTriAngle);

	unordered_map<size_t, Eigen::Vector3d> projectionCenters;
	for (size_t point3DID : points3DID)
	{
		if (points3D.find(point3DID) == points3D.end()) continue;
		const CPoint3D& point3D = points3D[point3DID];
		const vector<CTrackElement>& trackElements = point3D.GetTrack().GetAllElements();

		// 计算轨迹中所有影像姿态的两两组合的交会角, 仅当没有任何组合具有足够大的交会角时才删除点
		bool isKeepPoint = false;
		for (size_t i = 0; i < trackElements.size(); i++)
		{
			const size_t imageID1 = trackElements[i].imageID;
			Eigen::Vector3d projectionCenter1;
			if (projectionCenters.find(imageID1) == projectionCenters.end())
			{
				const CImage& image1 = database->GetImage(i);
				projectionCenter1 = image1.GetProjectionCenter(modelID);
				projectionCenters[imageID1] = projectionCenter1;
			}
			else
			{
				projectionCenter1 = projectionCenters[imageID1];
			}

			for (size_t j = 0; j < i; j++)
			{
				const size_t imageID2 = trackElements[j].imageID;
				const auto it = projectionCenters.find(imageID2);
				Check(it != projectionCenters.end());
				const Eigen::Vector3d projectionCenter2 = it->second;
				const double triAngle = CalculateTriangulationAngle(projectionCenter1, projectionCenter2, point3D.GetXYZ());
				if (triAngle >= minTriAngle_Rad)
				{
					isKeepPoint = true;
					break;
				}
			}
			if (isKeepPoint)
			{
				break;
			}
		}
		if (!isKeepPoint)
		{
			numFiltered++;
			DeletePoint3D(point3DID);
		}
	}
	return numFiltered;
}
size_t CModel::FilterPoints3DWithLargeReprojectionError(double maxReprojectionError, const unordered_set<size_t>& points3DID)
{
	Check(database);

	const double maxSquaredReprojectionError = maxReprojectionError * maxReprojectionError;
	size_t numFiltered = 0;
	for (size_t point3DID : points3DID)
	{
		if (points3D.find(point3DID) == points3D.end()) continue;

		CPoint3D& point3D = points3D[point3DID];
		const vector<CTrackElement>& trackElements = point3D.GetTrack().GetAllElements();
		if (trackElements.size() < 2)
		{
			numFiltered += trackElements.size();
			DeletePoint3D(point3DID);
			continue;
		}

		double reprojectionErrorSum = 0;
		vector<CTrackElement> trackElementsToDelete;
		for (const CTrackElement& trackElement : trackElements)
		{
			const CImage& image = database->GetImage(trackElement.imageID);
			const CCamera& camera = database->GetCamera(image.GetCameraID());
			const CKeypoint& point2D = image.GetKeypoint(trackElement.point2DIndex);
			const double squaredReprojectionError = CalculateSquaredReprojectionError(Eigen::Vector2d(point2D.pt.x, point2D.pt.y), point3D.GetXYZ(), image.GetWorldToCamera(modelID), camera);
			if (squaredReprojectionError > maxSquaredReprojectionError)
			{
				trackElementsToDelete.push_back(trackElement);
			}
			else
			{
				reprojectionErrorSum += sqrt(squaredReprojectionError);
			}
		}
		if (trackElementsToDelete.size() + 1 >= point3D.GetTrack().GetTrackLength())
		{
			numFiltered += point3D.GetTrack().GetTrackLength();
			DeletePoint3D(point3DID);
		}
		else
		{
			numFiltered += trackElementsToDelete.size();
			for (const CTrackElement& trackElement : trackElementsToDelete)
			{
				DeleteObservation(trackElement.imageID, trackElement.point2DIndex);
			}
			point3D.SetError(reprojectionErrorSum / point3D.GetTrack().GetTrackLength());
		}
	}
	return numFiltered;
}
tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> CModel::ComputeBoundsAndCentroid(double p0, double p1, bool isUseImages) const
{
	Check(p0 >= 0 && p0 <= 1 && p1 >= 0 && p1 <= 1 && p0 <= p1);
	const size_t numElements = (isUseImages ? regImageIDs.size() : points3D.size());
	if (numElements == 0)
	{
		return make_tuple(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0));
	}
	vector<float> coordsX;
	vector<float> coordsY;
	vector<float> coordsZ;
	if (isUseImages)
	{
		coordsX.reserve(regImageIDs.size());
		coordsY.reserve(regImageIDs.size());
		coordsZ.reserve(regImageIDs.size());
		for (size_t imageID : regImageIDs)
		{
			const Eigen::Vector3d projectionCenter = database->GetImage(imageID).GetProjectionCenter(modelID);
			coordsX.push_back(projectionCenter(0));
			coordsY.push_back(projectionCenter(1));
			coordsZ.push_back(projectionCenter(2));
		}
	}
	else
	{
		coordsX.reserve(points3D.size());
		coordsY.reserve(points3D.size());
		coordsZ.reserve(points3D.size());
		for (const auto& point3D : points3D) 
		{
			coordsX.push_back(point3D.second.GetX());
			coordsY.push_back(point3D.second.GetY());
			coordsZ.push_back(point3D.second.GetZ());
		}
	}

	sort(coordsX.begin(), coordsX.end());
	sort(coordsY.begin(), coordsY.end());
	sort(coordsZ.begin(), coordsZ.end());

	const size_t P0 = static_cast<size_t>((coordsX.size() > 3) ? p0 * (coordsX.size() - 1) : 0);
	const size_t P1 = static_cast<size_t>((coordsX.size() > 3) ? p1 * (coordsX.size() - 1) : coordsX.size() - 1);

	const Eigen::Vector3d bboxMin(coordsX[P0], coordsY[P0], coordsZ[P0]);
	const Eigen::Vector3d bboxMax(coordsX[P1], coordsY[P1], coordsZ[P1]);

	Eigen::Vector3d meanCoord(0, 0, 0);
	for (size_t i = P0; i <= P1; i++)
	{

		meanCoord(0) += coordsX[i];
		meanCoord(1) += coordsY[i];
		meanCoord(2) += coordsZ[i];
	}
	meanCoord /= P1 - P0 + 1;

	return make_tuple(bboxMin, bboxMax, meanCoord);
}





