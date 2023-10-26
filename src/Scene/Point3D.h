#pragma once
#include <opencv2/opencv.hpp>

#include "Track.h"
#include "../Base/Base.h"



class CPoint3D final
{
public:
	CPoint3D() noexcept
	{
		XYZ = Eigen::Vector3d(0, 0, 0);
		color = Eigen::Vector3ub(0, 0, 0);
		error = -1;
	}
	CPoint3D(double X, double Y, double Z) noexcept
	{
		XYZ = Eigen::Vector3d(X, Y, Z);
		color = Eigen::Vector3ub(0, 0, 0);
		error = -1;
	}
	CPoint3D(const Eigen::Vector3d& XYZ) noexcept
	{
		this->XYZ = XYZ;
		color = Eigen::Vector3ub(0, 0, 0);
		error = -1;
	}
	
	inline const Eigen::Vector3d& GetXYZ() const noexcept
	{
		return XYZ;
	}
	inline Eigen::Vector3d& GetXYZ() noexcept
	{
		return XYZ;
	}
	inline void GetXYZ(double& X, double& Y, double& Z) const noexcept
	{
		X = XYZ(0);
		Y = XYZ(1);
		Z = XYZ(2);
	}
	inline void SetXYZ(const Eigen::Vector3d& XYZ) noexcept
	{
		this->XYZ = XYZ;
	}
	inline void SetXYZ(double X, double Y, double Z) noexcept
	{
		XYZ(0) = X;
		XYZ(1) = Y;
		XYZ(2) = Z;
	}
	inline const double GetX() const noexcept
	{
		return XYZ(0);
	}
	inline const double GetY() const noexcept
	{
		return XYZ(1);
	}
	inline const double GetZ() const noexcept
	{
		return XYZ(2);
	}

	// RGB
	inline const Eigen::Vector3ub& GetColor() const noexcept
	{
		return color;
	}
	inline void GetColor(uchar& R, uchar& G, uchar& B) const noexcept
	{
		R = color[0];
		G = color[1];
		B = color[2];
	}
	inline void SetColor(const Eigen::Vector3ub& color_RGB) noexcept
	{
		color = color_RGB;
	}
	inline void SetColor(uchar R, uchar G, uchar B) noexcept
	{
		color[0] = R;
		color[1] = G;
		color[2] = B;
	}
	inline const uchar GetColor_R() const noexcept
	{
		return color[0];
	}
	inline void SetColor_R(uchar R) noexcept
	{
		color[0] = R;
	}
	inline const uchar GetColor_G() const noexcept
	{
		return color[1];
	}
	inline void SetColor_G(uchar G) noexcept
	{
		color[1] = G;
	}
	inline const uchar GetColor_B() const noexcept
	{
		return color[2];
	}
	inline void SetColor_B(uchar B) noexcept
	{
		color[2] = B;
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
	Eigen::Vector3d XYZ;
	Eigen::Vector3ub color;
	double error;
	CTrack track;
};