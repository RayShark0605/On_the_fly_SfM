#pragma once
#define EIGEN_USE_MKL_ALL
#include "Rigid3D.h"
#include "Sim3D.h"
#include "../Base/Math.h"

// 通过将给定矩阵的奇异值设为1, 计算具有最接近的Frobenius范数的最接近的旋转矩阵
Eigen::Matrix3d ComputeClosestRotationMatrix(const Eigen::Matrix3d& matrix);

// 将投影矩阵分解为内参矩阵, 旋转矩阵和平移向量. 如果分解失败, 返回false
bool DecomposeProjectionMatrix(const Eigen::Matrix3x4d& projectionMatrix, Eigen::Matrix3d& K, Eigen::Matrix3d& R, Eigen::Vector3d& T);

// 从一个向量构造反对称的叉积矩阵
Eigen::Matrix3d CrossProductMatrix(const Eigen::Vector3d& vector);

// 三维旋转矩阵转换成欧拉角. 使用的约定是: R = Rx * Ry * Rz, 并且使用右手系. 返回的rx, ry, rz为弧度制
void RotationMatrixToEulerAngles(const Eigen::Matrix3d& R, double& rx, double& ry, double& rz);

// 欧拉角转换成三维旋转矩阵. 使用的约定是: R = Rx * Ry * Rz, 并且使用右手系. 输入的rx, ry, rz必须为弧度制
Eigen::Matrix3d EulerAnglesToRotationMatrix(double rx, double ry, double rz);

// 计算多个四元数的加权平均. 参考文献: Markley, F. Landis, et al. "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30.4 (2007): 1193-1197.
Eigen::Quaterniond AverageQuaternions(const std::vector<Eigen::Quaterniond>& quats, const std::vector<double>& weights);

// 线性插值相机位姿
CRigid3D InterpolateCameraPoses(const CRigid3D& worldToCamera1, const CRigid3D& worldToCamera2, double t);

// 计算3D点在给定相机坐标系下的深度(Z坐标)
double CalculateDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D);

// 执行正景深约束测试, 即确定哪些三角化的对应点位于两个相机的前方. 第一个相机具有投影矩阵P1=[I | 0], 第二个相机具有投影矩阵P2=[R | t].
// @param R            第二个投影矩阵的3×3旋转矩阵
// @param t            第二个投影矩阵的3×1平移向量
// @param points1      影像1上的同名点.
// @param points2      影像2上的同名点.
// @param points3D     位于两个相机前方的3D点.
bool CheckCheirality(const Eigen::Matrix3d& R, const Eigen::Vector3d& t, const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, std::vector<Eigen::Vector3d>& points3D);

// 根据相机在旧世界坐标系下的位姿和新世界坐标系相对于旧世界坐标系的Sim3D变换, 计算相机在新世界坐标系下的位姿
CRigid3D TransformCameraWorld(const CSim3D& oldWorldToNewWorld, const CRigid3D& oldWorldToCamera);

// 使用2D-3D对应关系估计绝对姿态(可选地估计焦距)
// @param options         选项
// @param points2D        2D点
// @param points3D        3D点
// @param worldToCamera   估计的绝对相机姿态
// @param camera          需要估计姿态的相机. 在原地修改以存储估计的焦距
// @param numInliers      RANSAC中的内点数
// @param inlierMask      2D-3D对应关系的内点掩码
bool EstimateAbsolutePose(const COptions& options, const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, size_t& numInliers, std::vector<char>& inlierMask);


// 从2D-3D对应关系中精化绝对位姿(可选地精化焦距)
//
// @param options            选项
// @param inlierMask         2D-3D对应关系的内点掩码
// @param points2D           对应的2D点
// @param points3D           对应的3D点
// @param worldToCamera      精化后的绝对相机位姿
// @param camera             需要估计姿态的相机. 在原地修改以存储估计的焦距
// @param worldToCameraCov   6x6协方差矩阵(可选), 包括旋转(表现为坐标轴角度, 在切空间中)和平移项.
bool RefineAbsolutePose(const COptions& options, const std::vector<char>& inlierMask, const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, CRigid3D& worldToCamera, CCamera& camera, Eigen::Matrix6d* worldToCameraCov = nullptr);
