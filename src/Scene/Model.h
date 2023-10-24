#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Point3D.h"
#include "Database.h"

class CModel final
{
public:
	CModel(size_t modelID, CDatabase* database);

	inline size_t GetModelID() const
	{
		return modelID;
	}
	inline void SetModelID(size_t modelID)
	{
		this->modelID = modelID;
	}

	inline bool IsImageRegistered(size_t imageID) const
	{
		if (regImageIDs.find(imageID) == regImageIDs.end())
		{
			return false;
		}
		CHECK(database && database->IsImageExists(imageID));
		CHECK(database->GetImage(imageID).IsRegistered(modelID));
		return true;
	}
	inline size_t GetRegImagesNum() const noexcept
	{
		return regImageIDs.size();
	}
	
	inline size_t GetPoints3DNum() const noexcept
	{
		return points3D.size();
	}
	std::unordered_set<size_t> GetAllPoints3D() const;
	inline const CPoint3D& GetPoint3D(size_t point3DID) const
	{
		auto it = points3D.find(point3DID);
		CHECK(it != points3D.end());
		return it->second;
	}
	inline CPoint3D& GetPoint3D(size_t point3DID)
	{
		auto it = points3D.find(point3DID);
		CHECK(it != points3D.end());
		return it->second;
	}
	inline bool IsExistsPoint3D(size_t point3DID) const
	{
		return points3D.find(point3DID) != points3D.end();
	}
	size_t AddPoint3D(const Eigen::Vector3d& XYZ, const CTrack& track, const Eigen::Vector3ub& color);

	// Ϊ�Ѵ��ڵ�3D�������µĹ۲�
	void AddObservation(size_t point3DID, const CTrackElement& trackElement);

	// �ϲ�����3D�㲢�ҷ����µ�3D��ID. �ϲ���ԭ��: ��������3D���Track���Ƚ��м�Ȩƽ��
	size_t MergePoints3D(size_t point3D1ID, size_t point3D2ID);

	// ɾ��һ��3D��, �Լ��ڹ۲⵽��Ӱ�����������Ĺ���
	void DeletePoint3D(size_t point3DID);

	// ɾ��ĳ��Ӱ���ĳ��2D���������е�3D��֮��Ĺ۲�, ���������3D��ֻ�ܱ�2��Ӱ��۲⵽, ��ô��3D��Ҳ�ᱻɾ��
	void DeleteObservation(size_t imageID, size_t point2DID);

	// ע��Ӱ��
	void RegisterImage(size_t imageID);

	// ȡ��ע��Ӱ��
	void DeRegisterImage(size_t imageID);

	// ͨ�����ź�ƽ�ƶԳ������й�һ��, ������ƽ��֮������˻��Ŀ��ӻ�Ч��, ͬʱ����㷨����ֵ�ȶ���
	// ƽ��: ʹ��������Ļ��λ�õľ�ֵλ������ϵ��ԭ��
	// ����: ʹ����С�������������λ��extent��, p0��p1��ʾ�����ǵ�������ĵ���С�����ٷ�λ��
	void Normalize(double extent = 10.0, double p0 = 0.1, double p1 = 0.9, bool isUseImage = true);

	// ����3D�������
	Eigen::Vector3d ComputeCentroid(double p0 = 0.1, double p1 = 0.9) const;

	// ��������3D�����С�߽�
	std::pair<Eigen::Vector3d, Eigen::Vector3d> ComputeBoundingBox(double p0 = 0.0, double p1 = 1.0) const;

	// ������Ӱ���3D��Ӧ����ά�����Ա任
	void Transform(const CSim3D& oldWorldToNewWorld);

	// �ҳ��ڵ�ǰģ�ͺʹ���ģ���ж���ע���Ӱ��ID
	std::vector<size_t> FindCommonRegisteredImages(const CModel& otherModel) const;

	// ���˵���ͶӰ������, ����Ȼ��߽���ǹ�С��3D��
	size_t FilterPoints3D(double maxReprojectionError, double minTriAngle, const std::unordered_set<size_t>& pointsToBeFiltered);
	size_t FilterPoints3DInImages(double maxReprojectionError, double minTriAngle, const std::unordered_set<size_t>& imageIDs);
	size_t FilterAllPoints3D(double maxReprojectionError, double minTriAngle);

	// ���˵�����ȵĹ۲�
	size_t FilterObservationsWithNegativeDepth();

	// ���˵�û�й۲����������������Ӱ��
	std::vector<size_t> FilterImages(double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam);

	size_t ComputeNumObservations() const;
	double ComputeMeanTrackLength() const;
	double ComputeMeanObservationsPerRegImage() const;
	double ComputeMeanReprojectionError() const;

	void UpdatePoint3DErrors();

	bool ExtractColorsForImage(size_t imageID, const std::string& imagePath);

private:
	size_t modelID;
	CDatabase* const database;
	std::unordered_set<size_t> regImageIDs;
	std::unordered_map<size_t, CPoint3D> points3D;
	size_t nextPoint3DID;

	void SetObservationAsTriangulated(size_t imageID, size_t point2DID, bool isContinuedPoint3D);
	void ResetTriObservations(size_t imageID, size_t point2DID, bool isDeletedPoint3D);

	size_t FilterPoints3DWithSmallTriangulationAngle(double minTriAngle, const std::unordered_set<size_t>& points3DID);
	size_t FilterPoints3DWithLargeReprojectionError(double maxReprojectionError, const std::unordered_set<size_t>& points3DID);

	std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> ComputeBoundsAndCentroid(double p0, double p1, bool isUseImages) const;
};

















