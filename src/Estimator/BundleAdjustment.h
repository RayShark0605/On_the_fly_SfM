#pragma once
#include "../Base/Base.h"
#include "../Base/Options.h"
#include "../Scene/Model.h"
#include "../Scene/Database.h"

// �������ò���ƽ���Ӱ��, 3D������
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

	// �������ģ�͵Ĳв�����, �в���������ڹ۲���������
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

	// ������ӵ�Ӱ����������Ϊ���������. 
	// Ĭ�������, ���������Ӱ���������ǿɱ��. ע��, �ڵ�����Щ����֮ǰ, �����������Ӧ��ͼ��
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

	// ������ӵ�Ӱ�����̬����Ϊ���������, ��̬����ͶӰ�������ת��ƽ�Ʋ��ֶ����
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

	// ����λ�˵�ƽ�Ʋ���, ��˳�����̬������������[0, 1, 2]��Χ��, ���ұ�����Ψһ��. ע��, �ڵ�����Щ����֮ǰ, �����������Ӧ��ͼ��
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

// ����ceres-solver��ƽ����
class CBundleAdjuster
{
public:
	CBundleAdjuster(const COptions& options, const CBundleAdjustmentConfig& config, CDatabase* const database);

	bool Solve(CModel& model);

	// ��ȡ��һ��Solve֮��Ľ���ժҪ
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





















