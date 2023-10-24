#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <array>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "../Base/Base.h"
#include "Rigid3D.h"
#include "Pose.h"

// �����ʾ���E�ֽ�Ϊ�������ܵ���ת����(R1��R2)�Լ��������ܵ�ƽ������(t��-t)
void DecomposeEssentialMatrix(const Eigen::Matrix3d& E, Eigen::Matrix3d& R1, Eigen::Matrix3d& R2, Eigen::Vector3d& t);

// �Ӹ����ı��ʾ���E�зֽ�����ܵ�λ��. ��һ��Ӱ���λ�˼���ΪP=[I | 0]
void EssentialMatrixToPose(const Eigen::Matrix3d& E, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, Eigen::Matrix3d& R, Eigen::Vector3d& t, std::vector<Eigen::Vector3d>& points3D);

// ����������λ����ɱ��ʾ���E. �����һ��Ӱ���λ�˵�ͶӰ����P=[I | 0], ���ҵڶ���Ӱ���λ���Ǵ���������ϵ���������ϵ�ı任
Eigen::Matrix3d PoseToEssentialMatrix(const CRigid3D& camera1ToCamera2);

// Ѱ��"����Ӱ���", ������õ�ʽ: optimal_point1 ^ T * E * optimal_point2 = 0, EΪ���ʾ�����������
// �ο�����: Lindstrom, P., "Triangulation made easy", Computer Vision and Pattern Recognition (CVPR), 2010 IEEE Conference on , vol., no., pp.1554,1561, 13-18 June 2010
void FindOptimalImageObservations(const Eigen::Matrix3d& E, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2, Eigen::Vector2d& optimalPoint1, Eigen::Vector2d& optimalPoint2);

// ���㼫���λ��(����������ʾ). isLeftImageΪtrue��������Ӱ��ļ���, Ϊfalse��������Ӱ��ļ���
Eigen::Vector3d EssentialMatrixToEpipole(const Eigen::Matrix3d& E, bool isLeftImage);

// ���ʾ������. �������A�����B��ת��, �������B�����A��ת��
Eigen::Matrix3d InvertEssentialMatrix(const Eigen::Matrix3d& E);