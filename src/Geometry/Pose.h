#pragma once
#define EIGEN_USE_MKL_ALL
#include "Rigid3D.h"
#include "Sim3D.h"
#include "../Base/Math.h"

// ͨ�����������������ֵ��Ϊ1, ���������ӽ���Frobenius��������ӽ�����ת����
Eigen::Matrix3d ComputeClosestRotationMatrix(const Eigen::Matrix3d& matrix);

// ��ͶӰ����ֽ�Ϊ�ڲξ���, ��ת�����ƽ������. ����ֽ�ʧ��, ����false
bool DecomposeProjectionMatrix(const Eigen::Matrix3x4d& projectionMatrix, Eigen::Matrix3d& K, Eigen::Matrix3d& R, Eigen::Vector3d& T);

// ��һ���������췴�ԳƵĲ������
Eigen::Matrix3d CrossProductMatrix(const Eigen::Vector3d& vector);

// ��ά��ת����ת����ŷ����. ʹ�õ�Լ����: R = Rx * Ry * Rz, ����ʹ������ϵ. ���ص�rx, ry, rzΪ������
void RotationMatrixToEulerAngles(const Eigen::Matrix3d& R, double& rx, double& ry, double& rz);

// ŷ����ת������ά��ת����. ʹ�õ�Լ����: R = Rx * Ry * Rz, ����ʹ������ϵ. �����rx, ry, rz����Ϊ������
Eigen::Matrix3d EulerAnglesToRotationMatrix(double rx, double ry, double rz);

// ��������Ԫ���ļ�Ȩƽ��. �ο�����: Markley, F. Landis, et al. "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30.4 (2007): 1193-1197.
Eigen::Quaterniond AverageQuaternions(const std::vector<Eigen::Quaterniond>& quats, const std::vector<double>& weights);

// ���Բ�ֵ���λ��
CRigid3D InterpolateCameraPoses(const CRigid3D& worldToCamera1, const CRigid3D& worldToCamera2, double t);

// ����3D���ڸ����������ϵ�µ����(Z����)
double CalculateDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D);

// ִ��������Լ������, ��ȷ����Щ���ǻ��Ķ�Ӧ��λ�����������ǰ��. ��һ���������ͶӰ����P1=[I | 0], �ڶ����������ͶӰ����P2=[R | t].
// @param R            �ڶ���ͶӰ�����3��3��ת����
// @param t            �ڶ���ͶӰ�����3��1ƽ������
// @param points1      Ӱ��1�ϵ�ͬ����.
// @param points2      Ӱ��2�ϵ�ͬ����.
// @param points3D     λ���������ǰ����3D��.
bool CheckCheirality(const Eigen::Matrix3d& R, const Eigen::Vector3d& t, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, std::vector<Eigen::Vector3d>& points3D);

// ��������ھ���������ϵ�µ�λ�˺�����������ϵ����ھ���������ϵ��Sim3D�任, �������������������ϵ�µ�λ��
CRigid3D TransformCameraWorld(const CSim3D& oldWorldToNewWorld, const CRigid3D& oldWorldToCamera);

// ʹ��2D-3D��Ӧ��ϵ���ƾ�����̬(��ѡ�ع��ƽ���)
// @param options         ѡ��
// @param points2D        2D��
// @param points3D        3D��
// @param worldToCamera   ���Ƶľ��������̬
// @param camera          ��Ҫ������̬�����. ��ԭ���޸��Դ洢���ƵĽ���
// @param numInliers      RANSAC�е��ڵ���
// @param inlierMask      2D-3D��Ӧ��ϵ���ڵ�����
bool EstimateAbsolutePose(const COptions& options, const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, size_t& numInliers, std::vector<char>& inlierMask);


// ��2D-3D��Ӧ��ϵ�о�������λ��(��ѡ�ؾ�������)
//
// @param options            ѡ��
// @param inlierMask         2D-3D��Ӧ��ϵ���ڵ�����
// @param points2D           ��Ӧ��2D��
// @param points3D           ��Ӧ��3D��
// @param worldToCamera      ������ľ������λ��
// @param camera             ��Ҫ������̬�����. ��ԭ���޸��Դ洢���ƵĽ���
// @param worldToCameraCov   6x6Э�������(��ѡ), ������ת(����Ϊ������Ƕ�, ���пռ���)��ƽ����.
bool RefineAbsolutePose(const COptions& options, const std::vector<char>& inlierMask, const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, Eigen::Matrix6d* worldToCameraCov = nullptr);
