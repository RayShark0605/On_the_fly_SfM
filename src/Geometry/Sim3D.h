#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "../Base/Base.h"
#include "../Estimator/Estimator.h"

// 3D相似性变换, 包括7个自由度: 缩放(1), 旋转(3), 平移(3). 可以将a坐标系下的点集xa变换到b坐标系下的点集xb=scale*R*xa+t
struct CSim3D
{
	double scale = 1;
	Eigen::Quaterniond rotation = Eigen::Quaterniond::Identity();
	Eigen::Vector3d translation = Eigen::Vector3d::Zero();

	CSim3D() noexcept;
	CSim3D(double scale, const Eigen::Quaterniond& rotation, const Eigen::Vector3d& translation) noexcept;

	CSim3D Inverse() const; // 返回逆变换. 当前CSim3d不会变

	void SaveToFile(const std::string path) const;
	void ReadFromFile(const std::string path);

	Eigen::Matrix3x4d ToMatrix() const noexcept;

	// 估计从srcPoints到targetPoints的变换, 如何成功则返回true
	bool Estimate(const std::vector<Eigen::Vector3d>& srcPoints, const std::vector<Eigen::Vector3d>& targetPoints);
};

// 使用3D相似性变换去变换一个3D点: Xb = atob * Xa
// 注意: 如果要把Xa变换到Xb再变换到Xc, 应该这么写: Xc = btoc * (atob * Xa), 括号不可省略!
inline Eigen::Vector3d operator*(const CSim3D& atob, const Eigen::Vector3d& Xa)
{
	return atob.scale * (atob.rotation * Xa) + atob.translation;
}

// 连接多个3D相似性变换. 
// 计算a到d的Sim3d: atod = ctod * btoc * atob
inline CSim3D operator*(const CSim3D& btoc, const CSim3D& atob)
{
	CSim3D atoc;
	atoc.scale = btoc.scale * atob.scale;
	atoc.rotation = (btoc.rotation * atob.rotation).normalized();
	atoc.translation = btoc.translation + (btoc.scale * (btoc.rotation * atob.translation));
	return atoc;
}
