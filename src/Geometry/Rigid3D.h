#pragma once
#define EIGEN_USE_MKL_ALL
#include <Eigen/Geometry>
#include "../Base/Base.h"

// ��ά����任. ��6�����ɶ�(3����ת+3��ƽ��)
struct CRigid3D final
{
	Eigen::Quaterniond rotation = Eigen::Quaterniond::Identity(); // ��ת��Ԫ��
	Eigen::Vector3d translation = Eigen::Vector3d::Zero();        // ƽ������

	CRigid3D() = default;
	CRigid3D(const Eigen::Quaterniond& rotation, const Eigen::Vector3d& translation) noexcept
	{
		this->rotation = rotation;
		this->translation = translation;
	}

	// ����ת��Ԫ����ƽ������ת����3��4�еľ���
	inline Eigen::Matrix3x4d ToMatrix() const noexcept
	{
		Eigen::Matrix3x4d matrix;
		matrix.leftCols<3>() = rotation.toRotationMatrix();
		matrix.col(3) = translation;
		return matrix;
	}

	// ��ǰ����任����任("a��b" -> "b��a"��
	inline CRigid3D Inverse() const
	{
		CRigid3D inversed;
		inversed.rotation = rotation.inverse();
		inversed.translation = inversed.rotation * -translation;
		return inversed;
	}

	inline std::string WriteToString() const
	{
		std::ostringstream out;
		out << "{";
		out << rotation.w() << "," << rotation.x() << "," << rotation.y() << "," << rotation.z();
		out << "},{";
		out << translation.x() << "," << translation.y() << "," << translation.z();
		out << "}";
		return out.str();
	}

	inline void ReadFromString(const std::string& str)
	{
		std::istringstream in(str);
		std::vector<double> values;
		char ch;
		double num;

		in >> ch;
		CHECK(ch == '{');

		for (int i = 0; i < 4; i++)
		{
			in >> num;
			values.push_back(num);
			in >> ch;
		}
		rotation = Eigen::Quaterniond(values[0], values[1], values[2], values[3]);
		values.clear();

		in >> ch >> ch;
		CHECK(ch == '{');

		for (int i = 0; i < 3; i++)
		{
			in >> num;
			values.push_back(num);
			in >> ch;
		}
		translation = Eigen::Vector3d(values[0], values[1], values[2]);
	}

};

// ʹ�ø���任Rigid3Dȥ�任3D����x
// ����: ��һ�����ʽ��ʹ�ö���任ʱҪС��, ��ΪC++�еĳ˷�������Ǵ�������ֵ��, ����ʱ��Ҫ���������ָ����ȷ�Ķ���任��˳��
inline Eigen::Vector3d operator*(const CRigid3D& Rigid3D, const Eigen::Vector3d& x) noexcept
{
	return Rigid3D.rotation * x + Rigid3D.translation;
}

// ��AtoB�ĸ���任����Ӧ����һ������任BtoC, ���յõ�����任AtoC
inline CRigid3D operator*(const CRigid3D& BtoC, const CRigid3D& AtoB) noexcept
{
	CRigid3D result;
	result.rotation = (BtoC.rotation * AtoB.rotation).normalized();
	result.translation = BtoC.translation + (BtoC.rotation * AtoB.translation);
	return result;
}