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

using CConjugatePoints = std::unordered_map<size_t, size_t>; // CConjugatePoints[i]表示在影像i上的同名点(匹配)
using CObjectPoints = std::unordered_map<size_t, size_t>;    // CObjectPoints[i]表示在模型i上的3D点ID

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

	// 在重建之前构建三角化对应关系金字塔
	void Setup(const CCamera& camera, size_t point3DVisibilityPyramidLevel = 6);
	
	// 在重建之后还原三角化对应关系金字塔
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

	// 获取3D点(由2D点通过三角测量等方式得到的)的数量
	inline size_t GetNumPoints3D(size_t modelID) const
	{
		const auto it = numPoint3D.find(modelID);
		if (it == numPoint3D.end())
		{
			return 0;
		}
		return it->second;
	}

	// 获取能够看到3D点的观测数量(至少与另一影像中的一个三角测量点有对应关系的2D点的数量)
	inline size_t GetNumVisiblePoints3D(size_t modelID) const
	{
		const auto it = numVisiblePoint3D.find(modelID);
		if (it == numVisiblePoint3D.end())
		{
			return 0;
		}
		return it->second;
	}

	// 获取观测数("至少与另一张影像存在对应关系"的2D点的数量)
	inline size_t GetNumObservations() const noexcept
	{
		return numObservations;
	}
	inline void SetNumObservations(size_t numObservations) noexcept
	{
		this->numObservations = numObservations;
	}

	// 获取"该影像上的一个2D点与另一张影像上的一个2D点相关联"的响应数
	inline size_t GetNumCorrespondences() const noexcept
	{
		return numCorrespondences;
	}
	inline void SetNumCorrespondences(size_t numCorrespondences) noexcept
	{
		this->numCorrespondences = numCorrespondences;
	}

	// 获取三角测量观测值的得分. 与NumVisiblePoints3D不同, 得分还考虑了三角测量观测值的分布情况, 用来在增量式重建中选择下一个最佳影像.
	inline size_t GetPoint3DVisibilityScore() const noexcept
	{
		return point3DVisibilityPyramid.Score();
	}

	// 获取POS(世界坐标系到相机坐标系的刚体变换)
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

	// 获取POS的先验值(例如从EXIF陀螺仪标签得到)
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
	
	// 是否存在另一幅影像中的2D点与该2D点共同观测到对应的3D点
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

	// 设置modelID号模型的point3DIndex号3D点是point2DIndex号2D点的三角观测点(例如: 它是3D点跟踪轨迹的一部分)
	void SetPoint3DForPoint2D(size_t point2DIndex, size_t point3DIndex, size_t modelID);

	// 设置当前2D点没有三角化(例如: 它不是3D点跟踪轨迹的一部分)
	void ResetPoint3DForPoint2D(size_t point2DIndex, size_t modelID);

	// 是否存在某个2D点对应这个3D点
	bool HasPoint3D(size_t modelID, size_t point3DIndex) const;

	// 某个2D点是否存在对应的3D点
	bool IsPoint2DHasPoint3D(size_t modelID, size_t point2DIndex) const;

	// 指示存在另一幅影像有一个已三角化的2D点与该2D点存在对应关系. 这个函数只能在Setup之后调用
	void IncrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex);

	// 表示另一幅影像有一个不再被三角化的点, 并且与该2D点具有对应关系. 必须在之前对同一个point2DIndex调用过IncrementCorrespondenceHasPoint3D
	void DecrementCorrespondenceHasPoint3D(size_t modelID, size_t point2DIndex);

	// 删除该影像在modelID号模型上的所有数据
	void DeleteModelInfo(size_t modelID);

	// 获取影像在世界坐标系中的投影中心
	Eigen::Vector3d GetProjectionCenter(size_t modelID) const;

	// 获取影像的观察方向
	Eigen::Vector3d GetViewDirection(size_t modelID) const;

	// 表示当前影像的第thisImagePoint2DID号2D点, 与imageID号影像的第point2DID号2D点是同名点
	inline void AddCorrespondence(size_t thisImagePoint2DID, size_t imageID, size_t point2DID)
	{
		Check(thisImagePoint2DID < keypoints.size() && keypoints.size() == correspondences.size());
		if (correspondences[thisImagePoint2DID].first.empty()) // 当前影像的thisImagePoint2DID号2D点还没有任何一个响应
		{
			numObservations++;
		}
		correspondences[thisImagePoint2DID].first[imageID] = point2DID;
		numCorrespondences++;
	}

	// 获取响应数据
	const std::pair<CConjugatePoints, CObjectPoints>& GetCorrespondences(size_t point2DID) const;

	std::string WriteToString() const;

private:
	size_t ID;
	std::string name;
	size_t cameraID;
	std::atomic_size_t numObservations;      // "至少与另一张影像存在对应关系"的2D点的数量
	std::atomic_size_t numCorrespondences;   // 一个响应表示"该影像上的一个2D点与另一张影像上的一个2D点相关联"
	
	std::unordered_map<size_t, CRigid3D> worldToCamera;      // 影像的pos: 从世界坐标系到相机坐标系的刚体变换
	std::unordered_map<size_t, CRigid3D> worldToCameraPrior; // 影像的pos(先验估计值)

	std::unordered_set<size_t> registeredModels;          // 当前影像被哪些模型注册了
	std::unordered_map<size_t, size_t> numPoint3D;        // numPoint3D[0]=60表示当前影像在0号模型中可以观测到60个3D点
	std::unordered_map<size_t, size_t> numVisiblePoint3D; // numVisiblePoint3D[0]=60表示当前影像在0号模型中有60个2D点参与前方交会生成了3D点

	// correspondences[1].first[2]=3: 表示当前影像的第1号2D点与2号影像的3号2D点是同名点. correspondences[1].second[0]=34: 表示当前影像的第1号2D点在0号模型上的对应3D点ID是34
	std::vector<std::pair<CConjugatePoints, CObjectPoints>> correspondences;
	CKeypoints keypoints;
	CSIFTDescriptors descriptors;
	std::vector<std::unordered_map<size_t, size_t>> numCorrExistPoint3D;     // numCorrExistPoint3D[20][2]=10表示在第2号模型中, 当前影像的第20号2D点对应的3D点存在多少个2D相应点
	CVisibilityPyramid point3DVisibilityPyramid;            // 计算影像中三角化对应关系的分布
};

