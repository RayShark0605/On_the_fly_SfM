#pragma once
#include <vector>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "../Geometry/Rigid3D.h"
#include "Camera.h"

// ������ͶӰ���(3D����Ӱ���е�ͶӰ��Ӱ��۲��֮���ŷ����þ���), ���3D��λ���������, ��ô�����᷵��numeric_limits<double>::max()
double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera);
double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera);

// ����Ƕ����(��������ĵ�3D���γɵ�������۲�����֮��ĽǶ�)
double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera);
double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera);

// ʹ�ù�һ��Ӱ���������Ƕ����(��������ĵ�3D���γɵ�������۲�����֮��ĽǶ�)
double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera);
double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera);

// ���3D���Ƿ�����������Լ��, ���õ�λ�����ǰ��������ͼ��ƽ����
bool HasPointPositiveDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D);