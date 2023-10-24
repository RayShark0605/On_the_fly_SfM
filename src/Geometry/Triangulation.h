#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <Eigen/Dense>

#include "EssentialMatrix.h"
#include "../Base/Base.h"
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


