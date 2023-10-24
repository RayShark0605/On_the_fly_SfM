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


