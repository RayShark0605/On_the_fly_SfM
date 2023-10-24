#pragma once
#include <vector>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "../Geometry/Rigid3D.h"
#include "Camera.h"

// 计算重投影误差(3D点在影像中的投影与影像观测点之间的欧几里得距离), 如果3D点位于相机后面, 那么函数会返回numeric_limits<double>::max()
double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera);
double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera);

// 计算角度误差(从相机中心到3D点形成的视线与观测视线之间的角度)
double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera);
double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera);

// 使用归一化影像点来计算角度误差(从相机中心到3D点形成的视线与观测视线之间的角度)
double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera);
double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera);

// 检查3D点是否满足正景深约束, 即该点位于相机前方而非在图像平面上
bool HasPointPositiveDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D);