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

// 使用DLT(直接线性变换)做三角测量, 从两个不同视角的影像中的对应2D点来生成3D点
Eigen::Vector3d TriangulatePoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);

// 在相机位姿非常准确的前提下做更精确的三角测量. 前提条件: 相机位姿必须非常准确, 否则应该使用TriangulatePoint
// 参考论文: P. Lindstrom, "Triangulation Made Easy," IEEE Computer Vision and Pattern Recognition 2010, pp. 1554-1561, June 2010.
Eigen::Vector3d TriangulateOptimalPoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);

// 从多个不同视角的影像中的对应2D点来做三角测量(TriangulatePoint的最小二乘版本)
Eigen::Vector3d TriangulateMultiViewPoint(const std::vector<Eigen::Matrix3x4d>& worldToCameras, const std::vector<Eigen::Vector2d>& points2D);

// 从两张影像的多对2D点中三角测量多个3D点
std::vector<Eigen::Vector3d> TriangulatePoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2);

// 从两张影像的多对2D点中三角测量多个最优的3D点
std::vector<Eigen::Vector3d> TriangulateOptimalPoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2);

// 计算交会角, 返回结果是弧度制
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

	// 对影像观测点做前方交会, 包括: 新点的创建, 现有点的延续, 合并点(如果该影像连接了不同的轨迹).
	// 注意: 该影像必须已经在模型中注册, 并且已经有姿态
	size_t TriangulateImage(const COptions& options, size_t imageID);

	// 完成影像的前方交会. 尝试为还没有做前方交会的观测创建新的轨迹, 并尝试完成现有轨迹, 返回完成的观测数
	size_t CompleteImage(const COptions& options, size_t imageID);

	// 为3D点完成轨迹. 递归地尝试将观测添加到一个轨迹中, 该轨迹可能由于不准确的姿态等原因而未能前方交会, 返回完成的观测数
	size_t CompleteTracks(const COptions& options, const std::unordered_set<size_t>& points3DIDs);

	// 完成所有3D点的轨迹, 返回完成的观测数
	size_t CompleteAllTracks(const COptions& options);

	// 合并3D点的轨迹, 返回合并的观测数
	size_t MergeTracks(const COptions& options, const std::unordered_set<size_t>& points3DIDs);

	// 合并所有3D点的轨迹, 返回合并的观测数量
	size_t MergeAllTracks(const COptions& options);

	// 对欠重建的影像对重新做前方交会. 在重建发生漂移的情况下通常会执行这个. 如果triRatio小于retriangulateMinRatio则表明影像对是欠重建的
	size_t Retriangulate(const COptions& options);

	// 指出一个3D点已经被修改
	void AddModifiedPoint3D(size_t point3DID);
	
	// 获取自从上次ClearModifiedPoints3D之后已更新的3D点
	const std::unordered_set<size_t>& GetModifiedPoints3D();

	// 清除已更新的3D点
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
	// 清除非法相机参数和尝试合并的缓存
	void ClearCaches();

	// 传递地寻找与其他影像的对应关系
	size_t Find(const COptions& options, size_t imageID, size_t point2DID, size_t transitivity, std::vector<CCorrData>& corrsData);
	
	// 尝试从给定的对应关系创建一个新的3D点
	size_t Create(const COptions& options, const std::vector<CCorrData>& corrsData);
	
	// 尝试用给定的对应关系延续3D点
	size_t Continue(const COptions& options, const CCorrData& refCorrData, const std::vector<CCorrData>& corrsData);
	
	// 尝试将3D点与其任何对应的3D点合并
	size_t Merge(const COptions& options, size_t point3DID);
	
	// 尝试传递性地完成一个3D点的轨迹
	size_t Complete(const COptions& options, size_t point3DID);
	
	// 检查相机是否有非法参数, 并缓存结果
	bool IsCameraBogusParams(const COptions& options, const CCamera& camera);

	// 判断该影像的第point2DID号2D点的观测是否为双视观测, 必须满足: 该2D点只有一个响应, 而且其响应也只有该2D点这一个响应
	bool IsTwoViewObservation(const CImage& image, size_t point2DID) const;


	void FindTransitiveCorrespondences(size_t imageID, size_t point2DID, size_t transitivity, CConjugatePoints& corrs);

};















