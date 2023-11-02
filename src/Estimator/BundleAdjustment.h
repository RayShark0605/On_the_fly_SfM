#pragma once
#include "../Base/Base.h"
#include "../Base/Options.h"
#include "../Scene/Model.h"
#include "../Scene/Database.h"

// 用于配置参与平差的影像, 3D点和相机
class CBundleAdjustmentConfig final
{
public:
	CBundleAdjustmentConfig();

	inline size_t GetNumImages() const noexcept
	{
		return imageIDs.size();
	}
	inline size_t GetNumPoints() const noexcept
	{
		return variablePoint3DIDs.size() + constantPoint3DIDs.size();
	}
	inline size_t GetNumConstantCameraIntrinsics() const noexcept
	{
		return constantCameraIntrinsics.size();
	}
	inline size_t GetNumConstantCameraPoses() const noexcept
	{
		return constantCameraPoses.size();
	}
	inline size_t GetNumConstantCameraPositions() const noexcept
	{
		return constantCameraPositions.size();
	}
	inline size_t GetNumVariablePoints() const noexcept
	{
		return variablePoint3DIDs.size();
	}
	inline size_t GetNumConstantPoints() const noexcept
	{
		return constantPoint3DIDs.size();
	}

	// 计算给定模型的残差数量, 残差的数量等于观测数的两倍
	size_t GetNumResiduals(const CModel& model, const CDatabase* const database) const;

	inline void AddImage(size_t imageID)
	{
		imageIDs.insert(imageID);
	}
	inline bool IsExistImage(size_t imageID) const
	{
		return imageIDs.find(imageID) != imageIDs.end();
	}
	inline void RemoveImage(size_t imageID)
	{
		const auto it = imageIDs.find(imageID);
		Check(it != imageIDs.end());
		imageIDs.erase(it);
	}

	// 将已添加的影像的相机设置为常量或变量. 
	// 默认情况下, 所有已添加影像的相机都是可变的. 注意, 在调用这些方法之前, 必须先添加相应的图像
	inline void SetConstantCameraIntrinsics(size_t cameraID)
	{
		constantCameraIntrinsics.insert(cameraID);
	}
	inline void SetVariableCameraIntrinsics(size_t cameraID)
	{
		const auto it = constantCameraIntrinsics.find(cameraID);
		Check(it != constantCameraIntrinsics.end());
		constantCameraIntrinsics.erase(it);
	}
	inline bool IsConstantCameraIntrinsics(size_t cameraID) const
	{
		return constantCameraIntrinsics.find(cameraID) != constantCameraIntrinsics.end();
	}

	// 将已添加的影像的姿态设置为常量或变量, 姿态是由投影矩阵的旋转和平移部分定义的
	inline void SetConstantCameraPose(size_t imageID)
	{
		Check(imageIDs.find(imageID) != imageIDs.end());
		Check(constantCameraPositions.find(imageID) == constantCameraPositions.end());
		constantCameraPoses.insert(imageID);
	}
	inline void SetVariableCameraPose(size_t imageID)
	{
		const auto it = constantCameraPoses.find(imageID);
		Check(it != constantCameraPoses.end());
		constantCameraPoses.erase(it);
	}
	inline bool IsConstantCameraPose(size_t imageID) const
	{
		return constantCameraPoses.find(imageID) != constantCameraPoses.end();
	}

	// 设置位姿的平移部分, 因此常量姿态的索引可能在[0, 1, 2]范围内, 并且必须是唯一的. 注意, 在调用这些方法之前, 必须先添加相应的图像
	inline void SetConstantCameraPositions(size_t imageID, const std::vector<int>& IDs)
	{
		Check(!IDs.empty() && IDs.size() <= 3);
		Check(imageIDs.find(imageID) != imageIDs.end());
		Check(constantCameraPoses.find(imageID) == constantCameraPoses.end());
		Check(!IsVectorContainsDuplicateValues(IDs));

		constantCameraPositions.emplace(imageID, IDs);
	}
	inline void RemoveConstantCameraPositions(size_t imageID)
	{
		auto it = constantCameraPositions.find(imageID);
		Check(it != constantCameraPositions.end());
		constantCameraPositions.erase(it);
	}
	inline bool IsConstantCameraPositions(size_t imageID) const
	{
		return constantCameraPositions.find(imageID) != constantCameraPositions.end();
	}

	inline void AddVariablePoint(size_t point3DID)
	{
		Check(constantPoint3DIDs.find(point3DID) == constantPoint3DIDs.end());
		variablePoint3DIDs.insert(point3DID);
	}
	inline void AddConstantPoint(size_t point3DID)
	{
		Check(variablePoint3DIDs.find(point3DID) == variablePoint3DIDs.end());
		constantPoint3DIDs.insert(point3DID);
	}
	inline bool IsExistPoint(size_t point3DID) const
	{
		return variablePoint3DIDs.find(point3DID) != variablePoint3DIDs.end() || constantPoint3DIDs.find(point3DID) != constantPoint3DIDs.end();
	}
	inline bool IsExistVariablePoint(size_t point3DID) const
	{
		return variablePoint3DIDs.find(point3DID) != variablePoint3DIDs.end();
	}
	inline bool IsExistConstantPoint(size_t point3DID) const
	{
		return constantPoint3DIDs.find(point3DID) != constantPoint3DIDs.end();
	}
	inline void RemoveVariablePoint(size_t point3DID)
	{
		const auto it = variablePoint3DIDs.find(point3DID);
		Check(it != variablePoint3DIDs.end());
		variablePoint3DIDs.erase(it);
	}
	inline void RemoveConstantPoint(size_t point3DID)
	{
		const auto it = constantPoint3DIDs.find(point3DID);
		Check(it != constantPoint3DIDs.end());
		constantPoint3DIDs.erase(it);
	}

	inline const std::unordered_set<size_t>& GetAllImages() const
	{
		return imageIDs;
	}
	inline const std::unordered_set<size_t>& GetAllVariablePoints() const
	{
		return variablePoint3DIDs;
	}
	inline const std::unordered_set<size_t>& GetAllConstantPoints() const
	{
		return constantPoint3DIDs;
	}
	inline const std::vector<int>& GetConstantCameraPositions(size_t imageID) const
	{
		const auto it = constantCameraPositions.find(imageID);
		Check(it != constantCameraPositions.end());
		return it->second;
	}

private:
	std::unordered_set<size_t> constantCameraIntrinsics;
	std::unordered_set<size_t> imageIDs;
	std::unordered_set<size_t> variablePoint3DIDs;
	std::unordered_set<size_t> constantPoint3DIDs;
	std::unordered_set<size_t> constantCameraPoses;
	std::unordered_map<size_t, std::vector<int>> constantCameraPositions;
};

// 基于ceres-solver的平差器
class CBundleAdjuster
{
public:
	CBundleAdjuster(const COptions& options, const CBundleAdjustmentConfig& config, CDatabase* const database);

	bool Solve(CModel& model);

	// 获取上一次Solve之后的解算摘要
	inline const ceres::Solver::Summary& GetSummary() const noexcept
	{
		return summary;
	}



private:
	const COptions options;
	CDatabase* const database;
	CBundleAdjustmentConfig config;
	ceres::Problem* problem = nullptr;
	ceres::Solver::Summary summary;
	std::unordered_set<size_t> cameraIDs;
	std::unordered_map<size_t, size_t> point3DNumObservations;

	void Setup(CModel& model, ceres::LossFunction* lossFunction);
	void AddImage(size_t imageID, CModel& model, ceres::LossFunction* lossFunction);
	void AddPoint(size_t point3DID, CModel& model, ceres::LossFunction* lossFunction);
	void ParameterizeCameras(CModel& model);
	void ParameterizePoints(CModel& model);

};





















