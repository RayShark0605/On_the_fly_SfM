#pragma once
#include <string>
#ifndef Q_MOC_RUN
#if defined(emit)
#undef emit
#include <tbb/tbb.h>
#define emit // restore the macro definition of "emit", as it was defined in gtmetamacros.h
#else
#include <tbb/tbb.h>
#endif // defined(emit)
#endif // Q_MOC_RUN


#include "Point2D.h"
#include "VisibilityPyramid.h"
#include "Camera.h"
#include "../Geometry/Rigid3D.h"

using CConjugatePoints = std::unordered_map<size_t, size_t>; // CConjugatePoints[i]��ʾ��Ӱ��i�ϵ�ͬ����(ƥ��)
using CObjectPoints = std::unordered_map<size_t, size_t>;    // CObjectPoints[i]��ʾ��ģ��i�ϵ�3D��ID

class CImage final
{
public:
	CImage();

	bool operator==(const CImage& other) const
	{
		return name == other.name;
	}
	CImage(const CImage& other)
	{
		ID = other.ID;
		name = other.name;
		cameraID = other.cameraID;
		numObservations = other.numObservations.load();
		numCorrespondences = other.numCorrespondences.load();

		worldToCamera = other.worldToCamera;
		worldToCameraPrior = other.worldToCameraPrior;

		registeredModels = other.registeredModels;
		numPoint3D = other.numPoint3D;
		numVisiblePoint3D = other.numVisiblePoint3D;

		correspondences = other.correspondences;
		keypoints = other.keypoints;
		descriptors = other.descriptors;
		numCorrExistPoint3D = other.numCorrExistPoint3D;

		point3DVisibilityPyramid = other.point3DVisibilityPyramid;
	}
	CImage& operator=(const CImage& other)
	{
		if (this == &other)
		{
			return *this;
		}

		ID = other.ID;
		name = other.name;
		cameraID = other.cameraID;
		numObservations = other.numObservations.load();
		numCorrespondences = other.numCorrespondences.load();

		worldToCamera = other.worldToCamera;
		worldToCameraPrior = other.worldToCameraPrior;

		registeredModels = other.registeredModels;
		numPoint3D = other.numPoint3D;
		numVisiblePoint3D = other.numVisiblePoint3D;

		correspondences = other.correspondences;
		keypoints = other.keypoints;
		descriptors = other.descriptors;
		numCorrExistPoint3D = other.numCorrExistPoint3D;

		point3DVisibilityPyramid = other.point3DVisibilityPyramid;

		return *this;
	}

	// ���ؽ�֮ǰ�������ǻ���Ӧ��ϵ������
	void Setup(const CCamera& camera, size_t point3DVisibilityPyramidLevel = 6);
	
	// ���ؽ�֮��ԭ���ǻ���Ӧ��ϵ������
	void TearDown() noexcept;

	inline size_t GetImageID() const noexcept
	{
		return ID;
	}
	inline void SetImageID(size_t newID) noexcept
	{
		ID = newID;
	}

	inline std::string GetImageName() const noexcept
	{
		return name;
	}
	inline std::string& GetImageName()
	{
		return name;
	}
	inline void SetImageName(const std::string newName) noexcept
	{
		name = newName;
	}

	inline size_t GetCameraID() const noexcept
	{
		return cameraID;
	}
	inline void SetCameraID(size_t newID) noexcept
	{
		cameraID = newID;
	}
	inline bool HasCamera() const noexcept
	{
		return cameraID != std::numeric_limits<size_t>::max();
	}

	inline bool IsRegistered(size_t modelID) const
	{
		return registeredModels.find(modelID) != registeredModels.end();
	}
	inline void SetRegistered(size_t modelID, bool isRegistered)
	{
		const auto it = registeredModels.find(modelID);
		if (isRegistered)
		{
			Check(it == registeredModels.end());
			registeredModels.insert(modelID);
		}
		else
		{
			Check(it != registeredModels.end());
			registeredModels.erase(it);
		}
	}

	inline size_t GetNumPoints2D() const noexcept
	{
		return keypoints.size();
	}
	inline const CKeypoint& GetKeypoint(size_t point2DID) const
	{
		Check(point2DID < keypoints.size());
		return keypoints[point2DID];
	}
	inline CKeypoint& GetKeypoint(size_t point2DID)
	{
		Check(point2DID < keypoints.size());
		return keypoints[point2DID];
	}
	inline void SetPoints2D(const CKeypoints& keypoints)
	{
		this->keypoints = keypoints;
		correspondences.resize(keypoints.size(), { CConjugatePoints(), CObjectPoints() });
	}
	inline const CKeypoints& GetKeypoints() const
	{
		return keypoints;
	}
	inline CKeypoints& GetKeypoints()
	{
		return keypoints;
	}
	inline size_t GetNumDescriptors() const noexcept
	{
		return descriptors.rows();
	}
	inline void SetDescriptors(const CSIFTDescriptors& descriptors)
	{
		this->descriptors = descriptors;
	}
	inline const CSIFTDescriptors& GetDescriptors() const
	{
		return descriptors;
	}
	inline CSIFTDescriptors& GetDescriptors()
	{
		return descriptors;
	}
	
	inline size_t GetPoint3DID(size_t point2DID, size_t modelID) const
	{
		Check(point2DID < correspondences.size());
		const auto it = correspondences[point2DID].second.find(modelID);
		Check(it != correspondences[point2DID].second.end());
		return it->second;
	}

	// ��ȡ3D��(��2D��ͨ�����ǲ����ȷ�ʽ�õ���)������
	inline size_t GetNumPoints3D(size_t modelID) const
	{
		const auto it = numPoint3D.find(modelID);
		if (it == numPoint3D.end())
		{
			return 0;
		}
		return it->second;
	}

	// ��ȡ�ܹ�����3D��Ĺ۲�����(��������һӰ���е�һ�����ǲ������ж�Ӧ��ϵ��2D�������)
	inline size_t GetNumVisiblePoints3D(size_t modelID) const
	{
		const auto it = numVisiblePoint3D.find(modelID);
		if (it == numVisiblePoint3D.end())
		{
			return 0;
		}
		return it->second;
	}

	// ��ȡ�۲���("��������һ��Ӱ����ڶ�Ӧ��ϵ"��2D�������)
	inline size_t GetNumObservations() const noexcept
	{
		return numObservations;
	}
	inline void SetNumObservations(size_t numObservations) noexcept
	{
		this->numObservations = numObservations;
	}

	// ��ȡ"��Ӱ���ϵ�һ��2D������һ��Ӱ���ϵ�һ��2D�������"����Ӧ��
	inline size_t GetNumCorrespondences() const noexcept
	{
		return numCorrespondences;
	}
	inline void SetNumCorrespondences(size_t numCorrespondences) noexcept
	{
		this->numCorrespondences = numCorrespondences;
	}

	// ��ȡ���ǲ����۲�ֵ�ĵ÷�. ��NumVisiblePoints3D��ͬ, �÷ֻ����������ǲ����۲�ֵ�ķֲ����, ����������ʽ�ؽ���ѡ����һ�����Ӱ��.
	inline size_t GetPoint3DVisibilityScore() const noexcept
	{
		return point3DVisibilityPyramid.Score();
	}

	// ��ȡPOS(��������ϵ���������ϵ�ĸ���任)
	inline const CRigid3D& GetWorldToCamera(size_t modelID) const
	{
		Check(registeredModels.find(modelID) != registeredModels.end());
		const auto it = worldToCamera.find(modelID);
		if (it != worldToCamera.end())
		{
			return it->second;
		}
		return CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d::Zero());
	}
	inline CRigid3D& GetWorldToCamera(size_t modelID)
	{
		Check(registeredModels.find(modelID) != registeredModels.end());
		const auto it = worldToCamera.find(modelID);
		if (it != worldToCamera.end())
		{
			return it->second;
		}
		worldToCamera[modelID] = CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d::Zero());
		return worldToCamera[modelID];
	}

	// ��ȡPOS������ֵ(�����EXIF�����Ǳ�ǩ�õ�)
	inline const CRigid3D& GetWorldToCameraPrior(size_t modelID) const
	{
		Check(registeredModels.find(modelID) != registeredModels.end());
		const auto it = worldToCameraPrior.find(modelID);
		if (it != worldToCameraPrior.end())
		{
			return it->second;
		}
		constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
		return CRigid3D(Eigen::Quaterniond(NaN, NaN, NaN, NaN), Eigen::Vector3d(NaN, NaN, NaN));
	}
	inline CRigid3D& GetWorldToCameraPrior(size_t modelID)
	{
		Check(registeredModels.find(modelID) != registeredModels.end());
		const auto it = worldToCameraPrior.find(modelID);
		if (it != worldToCameraPrior.end())
		{
			return it->second;
		}
		constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
		worldToCameraPrior[modelID] = CRigid3D(Eigen::Quaterniond(NaN, NaN, NaN, NaN), Eigen::Vector3d(NaN, NaN, NaN));
		return worldToCameraPrior[modelID];
	}

	inline std::vector<std::unordered_map<size_t, size_t>> GetAllNumCorrExistPoint3D() const
	{
		return numCorrExistPoint3D;
	}
	inline std::vector<size_t> GetAllNumCorrExistPoint3D(size_t modelID) const
	{
		const size_t numPoint2D = numCorrExistPoint3D.size();
		std::vector<size_t> result(numPoint2D, 0);
		for (size_t i = 0; i < numPoint2D; i++)
		{
			const auto it = numCorrExistPoint3D[i].find(modelID);
			if (it != numCorrExistPoint3D[i].end())
			{
				result[i] = it->second;
			}
		}
		return result;
	}
	
	// �Ƿ������һ��Ӱ���е�2D�����2D�㹲ͬ�۲⵽��Ӧ��3D��
	inline bool IsPoint3DVisible(size_t modelID, size_t point2DIndex) const
	{
		Check(point2DIndex < numCorrExistPoint3D.size());
		const auto it = numCorrExistPoint3D[point2DIndex].find(modelID);
		if (it == numCorrExistPoint3D[point2DIndex].end())
		{
			return false;
		}
		return it->second > 0;
	}

	// ����modelID��ģ�͵�point3DIndex��3D����point2DIndex��2D������ǹ۲��(����: ����3D����ٹ켣��һ����)
	void SetPoint3DForPoint2D(size_t point2DIndex, size_t point3DIndex, size_t modelID);

	// ���õ�ǰ2D��û�����ǻ�(����: ������3D����ٹ켣��һ����)
	void ResetPoint3DForPoint2D(size_t point2DIndex, size_t modelID);

	// �Ƿ����ĳ��2D���Ӧ���3D��
	bool HasPoint3D(size_t modelID, size_t point3DIndex) const;

	// ĳ��2D���Ƿ���ڶ�Ӧ��3D��
	bool IsPoint2DHasPoint3D(size_t modelID, size_t point2DIndex) const;

	// ָʾ������һ��Ӱ����һ�������ǻ���2D�����2D����ڶ�Ӧ��ϵ. �������ֻ����Setup֮�����
	void IncrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex);

	// ��ʾ��һ��Ӱ����һ�����ٱ����ǻ��ĵ�, �������2D����ж�Ӧ��ϵ. ������֮ǰ��ͬһ��point2DIndex���ù�IncrementCorrespondenceHasPoint3D
	void DecrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex);

	// ɾ����Ӱ����modelID��ģ���ϵ���������
	void DeleteModelInfo(size_t modelID);

	// ��ȡӰ������������ϵ�е�ͶӰ����
	Eigen::Vector3d GetProjectionCenter(size_t modelID) const;

	// ��ȡӰ��Ĺ۲췽��
	Eigen::Vector3d GetViewDirection(size_t modelID) const;

	// ��ʾ��ǰӰ��ĵ�thisImagePoint2DID��2D��, ��imageID��Ӱ��ĵ�point2DID��2D����ͬ����
	inline void AddCorrespondence(size_t thisImagePoint2DID, size_t imageID, size_t point2DID)
	{
		Check(thisImagePoint2DID < keypoints.size() && keypoints.size() == correspondences.size());
		if (correspondences[thisImagePoint2DID].first.empty()) // ��ǰӰ���thisImagePoint2DID��2D�㻹û���κ�һ����Ӧ
		{
			numObservations++;
		}
		correspondences[thisImagePoint2DID].first[imageID] = point2DID;
		numCorrespondences++;
	}

	// ��ȡ��Ӧ����
	const std::pair<CConjugatePoints, CObjectPoints>& GetCorrespondences(size_t point2DID) const;

	std::string WriteToString() const;

private:
	size_t ID;
	std::string name;
	size_t cameraID;
	std::atomic_size_t numObservations;      // "��������һ��Ӱ����ڶ�Ӧ��ϵ"��2D�������
	std::atomic_size_t numCorrespondences;   // һ����Ӧ��ʾ"��Ӱ���ϵ�һ��2D������һ��Ӱ���ϵ�һ��2D�������"
	
	std::unordered_map<size_t, CRigid3D> worldToCamera;      // Ӱ���pos: ����������ϵ���������ϵ�ĸ���任
	std::unordered_map<size_t, CRigid3D> worldToCameraPrior; // Ӱ���pos(�������ֵ)

	std::unordered_set<size_t> registeredModels;          // ��ǰӰ����Щģ��ע����
	std::unordered_map<size_t, size_t> numPoint3D;        // numPoint3D[0]=60��ʾ��ǰӰ����0��ģ���п��Թ۲⵽60��3D��
	std::unordered_map<size_t, size_t> numVisiblePoint3D; // numVisiblePoint3D[0]=60��ʾ��ǰӰ����0��ģ������60��2D�����ǰ������������3D��

	// correspondences[1].first[2]=3: ��ʾ��ǰӰ��ĵ�1��2D����2��Ӱ���3��2D����ͬ����. correspondences[1].second[0]=34: ��ʾ��ǰӰ��ĵ�1��2D����0��ģ���ϵĶ�Ӧ3D��ID��34
	std::vector<std::pair<CConjugatePoints, CObjectPoints>> correspondences;
	CKeypoints keypoints;
	CSIFTDescriptors descriptors;
	std::vector<std::unordered_map<size_t, size_t>> numCorrExistPoint3D;     // numCorrExistPoint3D[20][2]=10��ʾ�ڵ�2��ģ����, ��ǰӰ��ĵ�20��2D���Ӧ��3D����ڶ��ٸ�2D��Ӧ��
	CVisibilityPyramid point3DVisibilityPyramid;            // ����Ӱ�������ǻ���Ӧ��ϵ�ķֲ�
};

