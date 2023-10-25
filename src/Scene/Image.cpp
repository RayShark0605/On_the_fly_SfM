#include "Image.h"

using namespace std;

CImage::CImage()
{
	ID = numeric_limits<size_t>::max();
	name = "";
	cameraID = numeric_limits<size_t>::max();
	registeredModels.clear();
	numPoint3D.clear();
	numVisiblePoint3D.clear();
	numObservations = 0;
	numCorrespondences = 0;

	numCorrExistPoint3D.clear();
	keypoints.clear();
	correspondences.clear();
}
void CImage::Setup(const CCamera& camera, size_t point3DVisibilityPyramidLevel)
{
	CHECK(camera.GetCameraID() == cameraID && point3DVisibilityPyramidLevel > 0 && keypoints.size() == correspondences.size());
	point3DVisibilityPyramid = CVisibilityPyramid(point3DVisibilityPyramidLevel, camera.GetWidth(), camera.GetHeight());
}
void CImage::TearDown() noexcept
{
	point3DVisibilityPyramid = CVisibilityPyramid();
}
void CImage::SetPoint3DForPoint2D(size_t point2DIndex, size_t point3DIndex, size_t modelID)
{
	CHECK(point3DIndex != numeric_limits<size_t>::max() && point2DIndex < keypoints.size() && keypoints.size() == correspondences.size() && keypoints.size() == descriptors.rows());
	if (correspondences[point2DIndex].second.find(modelID) == correspondences[point2DIndex].second.end())
	{
		numPoint3D[modelID]++;
	}
	correspondences[point2DIndex].second[modelID] = point3DIndex;
}
void CImage::ResetPoint3DForPoint2D(size_t point2DIndex, size_t modelID)
{
	CHECK(point2DIndex < keypoints.size() && keypoints.size() == correspondences.size() && keypoints.size() == descriptors.rows());
	if (correspondences[point2DIndex].second.find(modelID) != correspondences[point2DIndex].second.end())
	{
		CHECK(numPoint3D[modelID] > 0);
		correspondences[point2DIndex].second.erase(modelID);
		numPoint3D[modelID]--;
		if (numPoint3D[modelID] == 0)
		{
			numPoint3D.erase(modelID);
		}
	}
}
bool CImage::HasPoint3D(size_t modelID, size_t point3DIndex) const
{
	CHECK(point3DIndex != numeric_limits<size_t>::max());
	for (size_t i = 0; i < correspondences.size(); i++)
	{
		const auto it = correspondences[i].second.find(modelID);
		if (it == correspondences[i].second.end())
		{
			continue;
		}
		if (it->second == point3DIndex)
		{
			return true;
		}
	}
	return false;
}
bool CImage::IsPoint2DHasPoint3D(size_t modelID, size_t point2DIndex) const
{
	CHECK(point2DIndex < correspondences.size());
	const auto it = correspondences[point2DIndex].second.find(modelID);
	return it != correspondences[point2DIndex].second.end() && it->second != numeric_limits<size_t>::max();
}
void CImage::IncrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex)
{
	CHECK(point2DIndex < keypoints.size() && keypoints.size() == correspondences.size() && keypoints.size() == descriptors.rows());
	numCorrExistPoint3D[point2DIndex][modelID]++;
	if (numCorrExistPoint3D[point2DIndex][modelID] == 1)
	{
		numVisiblePoint3D[modelID]++;
	}
	CHECK(numVisiblePoint3D[modelID] <= numObservations);

	const CKeypoint& point2D = keypoints[point2DIndex];
	point3DVisibilityPyramid.SetPoint(point2D.pt.x, point2D.pt.y);
}
void CImage::DecrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex)
{
	CHECK(point2DIndex < keypoints.size() && keypoints.size() == correspondences.size() && keypoints.size() == descriptors.rows());
	CHECK(numCorrExistPoint3D[point2DIndex][modelID] > 0);

	numCorrExistPoint3D[point2DIndex][modelID]--;
	if (numCorrExistPoint3D[point2DIndex][modelID] == 0)
	{
		numCorrExistPoint3D[point2DIndex].erase(modelID);
		numVisiblePoint3D[modelID]--;
		if (numVisiblePoint3D[modelID] == 0)
		{
			numVisiblePoint3D.erase(modelID);
		}
	}
	CHECK(numVisiblePoint3D.find(modelID) == numVisiblePoint3D.end() || numVisiblePoint3D[modelID] <= numObservations);

	const CKeypoint& point2D = keypoints[point2DIndex];
	point3DVisibilityPyramid.ResetPoint(point2D.pt.x, point2D.pt.y);
}
void CImage::DeleteModelInfo(size_t modelID)
{
	const auto findRegisteredModels = registeredModels.find(modelID);
	if (findRegisteredModels != registeredModels.end())
	{
		registeredModels.erase(findRegisteredModels);
	}

	const auto findNumPoint3D = numPoint3D.find(modelID);
	if (findNumPoint3D != numPoint3D.end())
	{
		numPoint3D.erase(findNumPoint3D);
	}

	const auto findNumVisiblePoint3D = numVisiblePoint3D.find(modelID);
	if (findNumVisiblePoint3D != numVisiblePoint3D.end())
	{
		numVisiblePoint3D.erase(findNumVisiblePoint3D);
	}

	for (std::pair<CConjugatePoints, CObjectPoints>& pair : correspondences)
	{
		const auto it = pair.second.find(modelID);
		if (it != pair.second.end())
		{
			pair.second.erase(it);
		}
	}
}
Eigen::Vector3d CImage::GetProjectionCenter(size_t modelID) const 
{
	const CRigid3D& worldToCamera = GetWorldToCamera(modelID);
	return worldToCamera.rotation.inverse() * -worldToCamera.translation;
}
Eigen::Vector3d CImage::GetViewDirection(size_t modelID) const
{
	const CRigid3D& worldToCamera = GetWorldToCamera(modelID);
	return worldToCamera.rotation.toRotationMatrix().row(2);
}
const pair<CConjugatePoints, CObjectPoints>& CImage::GetCorrespondences(size_t point2DID) const
{
	CHECK(point2DID < keypoints.size() && keypoints.size() == correspondences.size() && keypoints.size() == descriptors.rows());
	return correspondences[point2DID];
}
string CImage::WriteToString() const
{
	return "";
}