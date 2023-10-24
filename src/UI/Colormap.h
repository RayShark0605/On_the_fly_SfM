#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>


#include <Eigen/Core>
#include "../Scene/Point3D.h"


class CPointColormap
{
public:
	CPointColormap();
	virtual ~CPointColormap() = default;
	virtual void Prepare() = 0;
	virtual Eigen::Vector4f ComputeColor(size_t point3DID, const CPoint3D& point3D) = 0;


private:
	float scale, min, max, range;
	float minQ, maxQ;
};


class CImageColormap
{
public:
	CImageColormap();
	virtual ~CImageColormap() = default;
	virtual void Prepare() = 0;
};