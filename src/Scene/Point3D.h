#pragma once
#include <opencv2/opencv.hpp>

#include "Track.h"
#include "../Base/Base.h"



class CPoint3D final : public cv::Point3d
{
public:
	CPoint3D() noexcept
	{
		x = 0;
		y = 0;
		z = 0;
		color_BGR = cv::Vec3b(0, 0, 0);
		error = -1;
	}
	CPoint3D(double X, double Y, double Z) noexcept
	{
		x = X;
		y = Y;
		z = Z;
		color_BGR = cv::Vec3b(0, 0, 0);
		error = -1;
	}
	CPoint3D(const Eigen::Vector3d& XYZ) noexcept
	{
		x = XYZ.x();
		y = XYZ.y();
		z = XYZ.z();
		color_BGR = cv::Vec3b(0, 0, 0);
		error = -1;
	}
	
	inline const Eigen::Vector3d GetXYZ() const noexcept
	{
		return Eigen::Vector3d(x, y, z);
	}
	inline void GetXYZ(double& X, double& Y, double& Z) const noexcept
	{
		X = x;
		Y = y;
		Z = z;
	}
	inline void SetXYZ(const Eigen::Vector3d& vector3d) noexcept
	{
		x = vector3d.x();
		y = vector3d.y();
		z = vector3d.z();
	}
	inline void SetXYZ(double X, double Y, double Z) noexcept
	{
		x = X;
		y = Y;
		z = Z;
	}
	inline const double GetX() const noexcept
	{
		return x;
	}
	inline const double GetY() const noexcept
	{
		return y;
	}
	inline const double GetZ() const noexcept
	{
		return z;
	}

	// RGB
	inline const Eigen::Vector3ub GetColor() const noexcept
	{
		return Eigen::Vector3ub(color_BGR[2], color_BGR[1], color_BGR[0]);
	}
	inline void GetColor(uchar& R, uchar& G, uchar& B) const noexcept
	{
		R = color_BGR[2];
		G = color_BGR[1];
		B = color_BGR[0];
	}
	inline void SetColor(const Eigen::Vector3ub& color_RGB) noexcept
	{
		color_BGR[0] = color_RGB.z();
		color_BGR[1] = color_RGB.y();
		color_BGR[2] = color_RGB.x();
	}
	inline void SetColor(uchar R, uchar G, uchar B) noexcept
	{
		color_BGR[0] = B;
		color_BGR[1] = G;
		color_BGR[2] = R;
	}
	inline const uchar GetColor_R() const noexcept
	{
		return color_BGR[2];
	}
	inline void SetColor_R(uchar R) noexcept
	{
		color_BGR[2] = R;
	}
	inline const uchar GetColor_G() const noexcept
	{
		return color_BGR[1];
	}
	inline void SetColor_G(uchar G) noexcept
	{
		color_BGR[1] = G;
	}
	inline const uchar GetColor_B() const noexcept
	{
		return color_BGR[0];
	}
	inline void SetColor_B(uchar B) noexcept
	{
		color_BGR[0] = B;
	}

	inline const double GetError() const noexcept
	{
		return error;
	}
	inline bool HasError() const noexcept
	{
		return std::abs(error + 1.0) > 1e-5;
	}
	inline void SetError(double error) noexcept
	{
		this->error = error;
	}

	inline const CTrack& GetTrack() const noexcept
	{
		return track;
	}
	inline CTrack& GetTrack() noexcept
	{
		return track;
	}
	inline void SetTrack(const CTrack& track) noexcept
	{
		this->track = track;
	}

private:
	cv::Vec3b color_BGR; // BGR
	double error;
	CTrack track;
};