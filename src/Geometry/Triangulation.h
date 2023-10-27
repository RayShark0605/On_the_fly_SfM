#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <Eigen/Dense>

#include "EssentialMatrix.h"
#include "../Base/Base.h"
#include "../Base/Options.h"
#include "../Scene/Point2D.h"
#include "../Scene/Point3D.h"
#include "../Estimator/RANSAC.h"
#include "../Estimator/Estimator.h"

// ʹ��DLT(ֱ�����Ա任)�����ǲ���, ��������ͬ�ӽǵ�Ӱ���еĶ�Ӧ2D��������3D��
Eigen::Vector3d TriangulatePoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);

// �����λ�˷ǳ�׼ȷ��ǰ����������ȷ�����ǲ���. ǰ������: ���λ�˱���ǳ�׼ȷ, ����Ӧ��ʹ��TriangulatePoint
// �ο�����: P. Lindstrom, "Triangulation Made Easy," IEEE Computer Vision and Pattern Recognition 2010, pp. 1554-1561, June 2010.
Eigen::Vector3d TriangulateOptimalPoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);

// �Ӷ����ͬ�ӽǵ�Ӱ���еĶ�Ӧ2D���������ǲ���(TriangulatePoint����С���˰汾)
Eigen::Vector3d TriangulateMultiViewPoint(const std::vector<Eigen::Matrix3x4d>& worldToCameras, const std::vector<Eigen::Vector2d>& points2D);

// ������Ӱ��Ķ��2D�������ǲ������3D��
std::vector<Eigen::Vector3d> TriangulatePoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2);

// ������Ӱ��Ķ��2D�������ǲ���������ŵ�3D��
std::vector<Eigen::Vector3d> TriangulateOptimalPoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2);

// ���㽻���, ���ؽ���ǻ�����
double CalculateTriangulationAngle(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const Eigen::Vector3d& point3D);
std::vector<double> CalculateTriangulationAngles(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const std::vector<Eigen::Vector3d>& points3D);

bool EstimateTriangulation(const CEstimateTriangulationOptions& options, const std::vector<CTriangulationPoint>& points, const std::vector<CTriangulationPose>& poses, std::vector<char>& inlierMask, Eigen::Vector3d& XYZ);

class CDatabase;
class CModel;
class CImage;
using CConjugatePoints = std::unordered_map<size_t, size_t>;
class CTriangulator final
{
public:
	CTriangulator(CDatabase* const database, CModel& model);

	// ��Ӱ��۲����ǰ������, ����: �µ�Ĵ���, ���е������, �ϲ���(�����Ӱ�������˲�ͬ�Ĺ켣).
	// ע��: ��Ӱ������Ѿ���ģ����ע��, �����Ѿ�����̬
	size_t TriangulateImage(const COptions& options, size_t imageID);

	// ���Ӱ���ǰ������. ����Ϊ��û����ǰ������Ĺ۲ⴴ���µĹ켣, ������������й켣, ������ɵĹ۲���
	size_t CompleteImage(const COptions& options, size_t imageID);

	// Ϊ3D����ɹ켣. �ݹ�س��Խ��۲���ӵ�һ���켣��, �ù켣�������ڲ�׼ȷ����̬��ԭ���δ��ǰ������, ������ɵĹ۲���
	size_t CompleteTracks(const COptions& options, const std::unordered_set<size_t>& points3DIDs);

	// �������3D��Ĺ켣, ������ɵĹ۲���
	size_t CompleteAllTracks(const COptions& options);

	// �ϲ�3D��Ĺ켣, ���غϲ��Ĺ۲���
	size_t MergeTracks(const COptions& options, const std::unordered_set<size_t>& points3DIDs);

	// �ϲ�����3D��Ĺ켣, ���غϲ��Ĺ۲�����
	size_t MergeAllTracks(const COptions& options);

	// ��Ƿ�ؽ���Ӱ���������ǰ������. ���ؽ�����Ư�Ƶ������ͨ����ִ�����. ���triRatioС��retriangulateMinRatio�����Ӱ�����Ƿ�ؽ���
	size_t Retriangulate(const COptions& options);

	// ָ��һ��3D���Ѿ����޸�
	void AddModifiedPoint3D(size_t point3DID);
	
	// ��ȡ�Դ��ϴ�ClearModifiedPoints3D֮���Ѹ��µ�3D��
	const std::unordered_set<size_t>& GetModifiedPoints3D();

	// ����Ѹ��µ�3D��
	void ClearModifiedPoints3D();

private:
	CModel& model;
	CDatabase* const database;
	const size_t modelID;
	std::unordered_map<size_t, bool> isCameraBogusParams;
	std::unordered_map<size_t, std::unordered_set<size_t>> mergeTrials;
	std::unordered_map<std::pair<size_t, size_t>, int, MatchPairHash, MatchPairEqual> retriangulateNumTrials;
	std::unordered_set<size_t> modifiedPoint3DIDs;

	struct CCorrData
	{
		size_t imageID;
		size_t point2DID;
		const CImage* image;
		const CCamera* camera;
		const CKeypoint* point2D;
	};
	// ����Ƿ���������ͳ��Ժϲ��Ļ���
	void ClearCaches();

	// ���ݵ�Ѱ��������Ӱ��Ķ�Ӧ��ϵ
	size_t Find(const COptions& options, size_t imageID, size_t point2DID, size_t transitivity, std::vector<CCorrData>& corrsData);
	
	// ���ԴӸ����Ķ�Ӧ��ϵ����һ���µ�3D��
	size_t Create(const COptions& options, const std::vector<CCorrData>& corrsData);
	
	// �����ø����Ķ�Ӧ��ϵ����3D��
	size_t Continue(const COptions& options, const CCorrData& refCorrData, const std::vector<CCorrData>& corrsData);
	
	// ���Խ�3D�������κζ�Ӧ��3D��ϲ�
	size_t Merge(const COptions& options, size_t point3DID);
	
	// ���Դ����Ե����һ��3D��Ĺ켣
	size_t Complete(const COptions& options, size_t point3DID);
	
	// �������Ƿ��зǷ�����, ��������
	bool IsCameraBogusParams(const COptions& options, const CCamera& camera);

	// �жϸ�Ӱ��ĵ�point2DID��2D��Ĺ۲��Ƿ�Ϊ˫�ӹ۲�, ��������: ��2D��ֻ��һ����Ӧ, ��������ӦҲֻ�и�2D����һ����Ӧ
	bool IsTwoViewObservation(const CImage& image, size_t point2DID) const;


	void FindTransitiveCorrespondences(size_t imageID, size_t point2DID, size_t transitivity, CConjugatePoints& corrs);

};















