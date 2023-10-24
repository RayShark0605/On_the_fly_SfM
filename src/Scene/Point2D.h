#pragma once
#include <limits>
#include <vector>
#include <limits>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>

#include "../Base/Base.h"

struct CPoint2D :public Eigen::Vector2d
{
	size_t point3DID = std::numeric_limits<size_t>::max();
	CPoint2D()
	{
		x() = 0;
		y() = 0;
		point3DID = std::numeric_limits<size_t>::max();
	}
	inline bool HasPoint3D() const noexcept
	{
		return (point3DID != std::numeric_limits<size_t>::max());
	}
	inline Eigen::Vector2d GetXY() const noexcept
	{
		return Eigen::Vector2d(x(), y());
	}
	inline std::string WriteToString() const
	{
		std::ostringstream out;
		out << "{";
		out << x() << "," << y();
		out << "},{";	
		out << point3DID;
		out << "}";
		return out.str();
	}
	inline void ReadFromString(const std::string& str)
	{
		std::istringstream in(str);
		char ch;
		double x_val, y_val;
		size_t id;

		in >> ch;
		CHECK(ch == '{');

		in >> x_val >> ch >> y_val;
		CHECK(ch == ',');

		in >> ch >> ch >> ch;
		CHECK(ch == '{');

		in >> id;

		x() = x_val;
		y() = y_val;
		point3DID = id;
	}

};

struct CKeypoint final :public cv::KeyPoint
{
	// �ؼ���ķ�����״. a11��a22����x���y���ϵ�����, a12��a21���Ƽ��к���ת
	/*
	* | a11 a12 |   | scaleX*cos(angle)  -scaleY*sin(angle+shear) |
	* |         | = |                                             |
	* | a21 a22 |	| scaleX*sin(angle)  scaleY*cos(angle+shear)  |
	*/
	float a11, a12, a21, a22;

	inline CKeypoint() noexcept
	{
		pt.x = 0;
		pt.y = 0;
		size = 1;  // scale�߶�Ĭ��Ϊ1
		angle = 0; // �Ƕ�Ĭ��Ϊ0

		a11 = 1;
		a12 = 0;
		a21 = 0;
		a22 = 1;
	}
	inline CKeypoint(float x, float y) noexcept
	{
		pt.x = x;
		pt.y = y;
		size = 1;  // scale�߶�Ĭ��Ϊ1
		angle = 0; // �Ƕ�Ĭ��Ϊ0

		a11 = 1;
		a12 = 0;
		a21 = 0;
		a22 = 1;
	}
	inline CKeypoint(float x, float y, float scale, float orientation)
	{
		CHECK(scale >= 0);
		pt.x = x;
		pt.y = y;
		size = scale;        // scale�߶�
		angle = orientation; // �Ƕ�

		const float scaleCosOrientation = scale * std::cos(orientation);
		const float scaleSinOrientation = scale * std::sin(orientation);
		a11 = scaleCosOrientation;
		a12 = -scaleSinOrientation;
		a21 = scaleSinOrientation;
		a22 = scaleCosOrientation;
	}
	inline CKeypoint(float x, float y, float a11, float a12, float a21, float a22) noexcept
	{
		pt.x = x;
		pt.y = y;
		this->a11 = a11;
		this->a12 = a12;
		this->a21 = a21;
		this->a22 = a22;

		size = GetScale();
		angle = GetOrientation();
	}

	// ��������, ����, ��ת, �������������ɹؼ���
	inline void FromShapeParams(float x, float y, float scaleX, float scaleY, float orientation, float shear)
	{
		CHECK(scaleX >= 0 && scaleY >= 0);
		pt.x = x;
		pt.y = y;
		a11 = scaleX * std::cos(orientation);
		a12 = -scaleY * std::sin(orientation + shear);
		a21 = scaleX * std::sin(orientation);
		a22 = scaleY * std::cos(orientation + shear);
		size = GetScale();
		angle = GetOrientation();
	}

	inline void Rescale(float scale)
	{
		CHECK(scale > 0);
		pt.x *= scale;
		pt.y *= scale;
		a11 *= scale;
		a12 *= scale;
		a21 *= scale;
		a22 *= scale;

		size = GetScale();
		angle = GetOrientation();
	}
	inline void Rescale(float scaleX, float scaleY)
	{
		CHECK(scaleX > 0 && scaleY > 0);
		pt.x *= scaleX;
		pt.y *= scaleY;
		a11 *= scaleX;
		a12 *= scaleY;
		a21 *= scaleX;
		a22 *= scaleY;

		size = GetScale();
		angle = GetOrientation();
	}
	inline float GetScaleX() const noexcept
	{
		return std::sqrtf(a11 * a11 + a21 * a21);
	}
	inline float GetScaleY() const noexcept
	{
		return std::sqrtf(a12 * a12 + a22 * a22);
	}
	inline float GetScale() const noexcept
	{
		const float scaleX = std::sqrtf(a11 * a11 + a21 * a21);
		const float scaleY = std::sqrtf(a12 * a12 + a22 * a22);
		return (scaleX + scaleY) / 2.0;
	}
	inline float GetOrientation() const noexcept
	{
		return std::atan2(a21, a11);
	}
	inline float GetShear() const noexcept
	{
		return std::atan2(-a12, a22) - std::atan2(a21, a11);
	}

	inline std::string WriteToString() const
	{
		std::ostringstream out;
		// ��� pt.x, pt.y
		out << "{";
		out << pt.x << "," << pt.y;
		out << "},{";
		// ��� size, angle
		out << size << "," << angle;
		out << "},{";
		// ��� a11, a12, a21, a22
		out << a11 << "," << a12 << "," << a21 << "," << a22;
		out << "}";
		return out.str();
	}
	inline void ReadFromString(const std::string& str)
	{
		std::istringstream in(str);
		char ch;
		float x_val, y_val, sz, ang, a11_val, a12_val, a21_val, a22_val;

		// ��ȡ��һ�� '{'
		in >> ch;
		CHECK(ch == '{');

		// ��ȡ pt.x, pt.y
		in >> x_val >> ch >> y_val;
		CHECK(ch == ',');

		// ��ȡ��һ�� '}'
		in >> ch;
		CHECK(ch == '}');

		// �������Ų���ȡ�ڶ��� '{'
		in >> ch >> ch;
		CHECK(ch == '{');

		// ��ȡ size, angle
		in >> sz >> ch >> ang;
		CHECK(ch == ',');

		// ��ȡ�ڶ��� '}'
		in >> ch;
		CHECK(ch == '}');

		// �������Ų���ȡ������ '{'
		in >> ch >> ch;
		CHECK(ch == '{');

		// ��ȡ a11, a12, a21, a22
		in >> a11_val >> ch >> a12_val >> ch >> a21_val >> ch >> a22_val;
		CHECK(ch == ',');

		// ��ȡ������ '}'
		in >> ch;
		CHECK(ch == '}');

		// ��ֵ
		pt.x = x_val;
		pt.y = y_val;
		size = sz;
		angle = ang;
		a11 = a11_val;
		a12 = a12_val;
		a21 = a21_val;
		a22 = a22_val;
	}
};
using CKeypoints = std::vector<CKeypoint>;
inline std::vector<cv::KeyPoint> CKeypointsToCVKeypoints(const CKeypoints& keypoints)
{
	std::vector<cv::KeyPoint> results(keypoints.size());
	for (size_t i = 0; i < keypoints.size(); i++)
	{
		results[i] = keypoints[i];
	}
	return results;
}


using CSIFTDescriptor = Eigen::Matrix<unsigned char, 1, 128, Eigen::RowMajor>;
using CSIFTDescriptors = Eigen::Matrix<unsigned char, Eigen::Dynamic, 128, Eigen::RowMajor>;
using CSIFTDescriptors_Float = Eigen::Matrix<float, Eigen::Dynamic, 128, Eigen::RowMajor>;

struct CSIFTMatch
{
	size_t point2DIndex1 = std::numeric_limits<size_t>::max();
	size_t point2DIndex2 = std::numeric_limits<size_t>::max();
	CSIFTMatch() {}
	CSIFTMatch(size_t point2DIndex1, size_t point2DIndex2)
	{
		this->point2DIndex1 = point2DIndex1;
		this->point2DIndex2 = point2DIndex2;
	}
	inline std::string WriteToString() const
	{
		std::ostringstream out;
		out << "{";
		out << point2DIndex1 << "," << point2DIndex2;
		out << "}";
		return out.str();
	}
	inline void ReadFromString(const std::string& str)
	{
		std::istringstream in(str);
		char ch;
		size_t idx1, idx2;

		in >> ch;
		CHECK(ch == '{');

		in >> idx1 >> ch >> idx2;
		CHECK(ch == ',');

		in >> ch;
		CHECK(ch == '}');

		point2DIndex1 = idx1;
		point2DIndex2 = idx2;
	}
};
using CSIFTMatches = std::vector<CSIFTMatch>;