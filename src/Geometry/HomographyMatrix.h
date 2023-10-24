#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <array>
#include <Eigen/Core>
#include <Eigen/Dense>

// ����Ӧ����H�ֽ�ɿ��ܵ���תR, ƽ��t��ƽ�淨����n. �ο�����: Malis, Ezio, and Manuel Vargas. "Deeper understanding of the homography decomposition for vision-based control." (2007): 90.
// ��һ��λ�˼���ΪP=[I | 0]. ע��: ������ص�R,t��n��sizeΪ4, ��˵����Ӧ������ƽ���յ���. ������ص�R,t��n��sizeΪ1, ��˵����Ӧ���Ǵ���ת��.
// 
// @param H          3x3��Ӧ����
// @param K1         ���һ���ڲξ���
// @param K2         ��������ڲξ���
// @param R          ���: ���п��ܵ�3��3��ת����
// @param t          ���: ���п��ܵ�ƽ������
// @param n          ���: ���п��ܵķ�����
void DecomposeHomographyMatrix(const Eigen::Matrix3d& H, const Eigen::Matrix3d& K1, const Eigen::Matrix3d& K2, std::vector<Eigen::Matrix3d>& R, std::vector<Eigen::Vector3d>& t, std::vector<Eigen::Vector3d>& n);

// �ӵ�Ӧ����H�ָ����п��ܵ�λ��(��һ��λ�˼���ΪP=[I | 0]). 
// 
// @param H          3x3��Ӧ����
// @param K1         ���һ���ڲξ���
// @param K2         ��������ڲξ���
// @param points1    Ӱ��һ�ϵĶ�Ӧ��
// @param points2    Ӱ����ϵĶ�Ӧ��
// @param R          ���: ���п��ܵ�3��3��ת����
// @param t          ���: ���п��ܵ�ƽ������
// @param n          ���: ���п��ܵķ�����
// @param points3D   ���: λ�����ǰ�������ǻ����3D��(������Ӧ�Բ��Ǵ���תʱ)
void HomographyMatrixToPose(const Eigen::Matrix3d& H, const Eigen::Matrix3d& K1, const Eigen::Matrix3d& K2, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, Eigen::Matrix3d& R, Eigen::Vector3d& t, Eigen::Vector3d& n, std::vector<Eigen::Vector3d>& points3D);


// �������λ�˹�����Ӧ����
//
// @param K1         ���һ���ڲξ���
// @param K2         ��������ڲξ���
// @param R          3��3��ת����
// @param t          ƽ������
// @param n          ������
// @param d          ��ƽ��Ĵ�ֱ����
//
// @return           3x3��Ӧ����
Eigen::Matrix3d PoseToHomographyMatrix(const Eigen::Matrix3d& K1, const Eigen::Matrix3d& K2, const Eigen::Matrix3d& R, const Eigen::Vector3d& t, const Eigen::Vector3d& n, double d);
