#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Sampler.h"
#include "SupportMeasurer.h"
#include "../Base/Base.h"


// 中心化和归一化图像上的点. 
// 会经过两步变换: 1. 中心化: 新的坐标系统的原点将位于图像点的质心. 2. 归一化: 使得点到新坐标系统原点的平均距离为sqrt(2)
//
// @param points            点坐标.
// @param normedPoints      输出: 转换后的点坐标.
// @param originToResult    输出: 3×3的转换矩阵.
void CenterAndNormalizeImagePoints(const std::vector<Eigen::Vector2d>& points, std::vector<Eigen::Vector2d>& normedPoints, Eigen::Matrix3d& originToResult);

// 计算一组对应点和给定的基础矩阵或本质矩阵的残差. 残差被定义为Sampson误差的平方
//
// @param points1     第一组对应点.
// @param points2     第二组对应点.
// @param E           3×3的基础矩阵或本质矩阵
// @param residuals   输出: 残差向量.
void ComputeSquaredSampsonError(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, std::vector<double>& residuals);

// 计算一组2D图像点和与之对应的3D点的平方重投影误差. 如果某个3D点在相机的后面, 那么它对应的平方重投影误差将被设置为std::numeric_limits<double>::max()
//
// @param points2D      归一化的2D图像点
// @param points3D      3D点
// @param proj_matrix   3×4的投影矩阵
// @param residuals     输出: 平方重投影误差向量
void ComputeSquaredReprojectionError(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& worldToCamera, std::vector<double>& residuals);
