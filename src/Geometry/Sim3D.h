#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "../Base/Base.h"
#include "../Estimator/Estimator.h"

// 3D�����Ա任, ����7�����ɶ�: ����(1), ��ת(3), ƽ��(3). ���Խ�a����ϵ�µĵ㼯xa�任��b����ϵ�µĵ㼯xb=scale*R*xa+t
struct CSim3D
{
	double scale = 1;
	Eigen::Quaterniond rotation = Eigen::Quaterniond::Identity();
	Eigen::Vector3d translation = Eigen::Vector3d::Zero();

	CSim3D() noexcept;
	CSim3D(double scale, const Eigen::Quaterniond& rotation, const Eigen::Vector3d& translation) noexcept;

	CSim3D Inverse() const; // ������任. ��ǰCSim3d�����

	void SaveToFile(const std::string path) const;
	void ReadFromFile(const std::string path);

	Eigen::Matrix3x4d ToMatrix() const noexcept;

	// ���ƴ�srcPoints��targetPoints�ı任, ��γɹ��򷵻�true
	bool Estimate(const std::vector<Eigen::Vector3d>& srcPoints, const std::vector<Eigen::Vector3d>& targetPoints);
};

// ʹ��3D�����Ա任ȥ�任һ��3D��: Xb = atob * Xa
// ע��: ���Ҫ��Xa�任��Xb�ٱ任��Xc, Ӧ����ôд: Xc = btoc * (atob * Xa), ���Ų���ʡ��!
inline Eigen::Vector3d operator*(const CSim3D& atob, const Eigen::Vector3d& Xa)
{
	return atob.scale * (atob.rotation * Xa) + atob.translation;
}

// ���Ӷ��3D�����Ա任. 
// ����a��d��Sim3d: atod = ctod * btoc * atob
inline CSim3D operator*(const CSim3D& btoc, const CSim3D& atob)
{
	CSim3D atoc;
	atoc.scale = btoc.scale * atob.scale;
	atoc.rotation = (btoc.rotation * atob.rotation).normalized();
	atoc.translation = btoc.translation + (btoc.scale * (btoc.rotation * atob.translation));
	return atoc;
}
