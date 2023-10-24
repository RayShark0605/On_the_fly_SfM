#pragma once
#define EIGEN_USE_MKL_ALL
#include <Eigen/Geometry>
#include "../Base/Base.h"

// 三维刚体变换. 有6个自由度(3个旋转+3个平移)
struct CRigid3D final
{
	Eigen::Quaterniond rotation = Eigen::Quaterniond::Identity(); // 旋转四元数
	Eigen::Vector3d translation = Eigen::Vector3d::Zero();        // 平移向量

	CRigid3D() = default;
	CRigid3D(const Eigen::Quaterniond& rotation, const Eigen::Vector3d& translation) noexcept
	{
		this->rotation = rotation;
		this->translation = translation;
	}

	// 将旋转四元数和平移向量转换成3行4列的矩阵
	inline Eigen::Matrix3x4d ToMatrix() const noexcept
	{
		Eigen::Matrix3x4d matrix;
		matrix.leftCols<3>() = rotation.toRotationMatrix();
		matrix.col(3) = translation;
		return matrix;
	}

	// 求当前刚体变换的逆变换("a到b" -> "b到a"）
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

// 使用刚体变换Rigid3D去变换3D向量x
// 警告: 在一个表达式中使用多个变换时要小心, 因为C++中的乘法运算符是从左到右求值的, 这种时候要添加括号来指定正确的多个变换的顺序
inline Eigen::Vector3d operator*(const CRigid3D& Rigid3D, const Eigen::Vector3d& x) noexcept
{
	return Rigid3D.rotation * x + Rigid3D.translation;
}

// 在AtoB的刚体变换上再应用另一个刚体变换BtoC, 最终得到刚体变换AtoC
inline CRigid3D operator*(const CRigid3D& BtoC, const CRigid3D& AtoB) noexcept
{
	CRigid3D result;
	result.rotation = (BtoC.rotation * AtoB.rotation).normalized();
	result.translation = BtoC.translation + (BtoC.rotation * AtoB.translation);
	return result;
}