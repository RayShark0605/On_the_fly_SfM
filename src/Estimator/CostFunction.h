#pragma once
#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "../Geometry/Rigid3D.h"
#include "../Scene/CameraModel.h"

template <typename T>
using EigenVector3Map = Eigen::Map<const Eigen::Matrix<T, 3, 1>>;
template <typename T>
using EigenQuaternionMap = Eigen::Map<const Eigen::Quaternion<T>>;

// 用于可变相机位姿, 标定参数和点参数的标准平差代价函数
class CReprojectionErrorCostFunction
{
public:
	explicit CReprojectionErrorCostFunction(const Eigen::Vector2d& point2D) : observedX(point2D(0)), observedY(point2D(1)){}
	static ceres::CostFunction* Create(const Eigen::Vector2d& point2D)
	{
		CReprojectionErrorCostFunction* reprojectionErrorCostFunction = new CReprojectionErrorCostFunction(point2D);
		ceres::CostFunction* costFunction = new ceres::AutoDiffCostFunction<CReprojectionErrorCostFunction, 2, 4, 3, 3, 4>(reprojectionErrorCostFunction);
		return costFunction;
	}
	template <typename T>
	bool operator()(const T* const worldToCameraRotation, const T* const worldToCameraTranslation, const T* const point3D, const T* const cameraParams, T* residuals) const
	{
		const Eigen::Matrix<T, 3, 1> point3DInCamera = EigenQuaternionMap<T>(worldToCameraRotation) * EigenVector3Map<T>(point3D) + EigenVector3Map<T>(worldToCameraTranslation);
		CSimpleRadialCameraModel simpleRadialCameraModel;
		simpleRadialCameraModel.CameraToImage(cameraParams, point3DInCamera[0], point3DInCamera[1], point3DInCamera[2], &residuals[0], &residuals[1]);
		residuals[0] -= T(observedX);
		residuals[1] -= T(observedY);
		return true;
	}

private:
	const double observedX;
	const double observedY;
};

// 用于可变相机标定参数和点参数, 固定相机位姿的平差代价函数
class CReprojectionErrorConstantPoseCostFunction
{
public:
	CReprojectionErrorConstantPoseCostFunction(const CRigid3D& worldToCamera, const Eigen::Vector2d& point2D) : worldToCamera(worldToCamera), observedX(point2D(0)), observedY(point2D(1)) {}
	static ceres::CostFunction* Create(const CRigid3D& worldToCamera, const Eigen::Vector2d& point2D)
	{
		CReprojectionErrorConstantPoseCostFunction* reprojectionErrorConstantPoseCostFunction = new CReprojectionErrorConstantPoseCostFunction(worldToCamera, point2D);
		ceres::CostFunction* costFunction = new ceres::AutoDiffCostFunction<CReprojectionErrorConstantPoseCostFunction, 2, 3, 4>(reprojectionErrorConstantPoseCostFunction);
		return costFunction;
	}
	template <typename T>
	bool operator()(const T* const point3D, const T* const cameraParams, T* residuals) const
	{
		const Eigen::Matrix<T, 3, 1> point3DInCamera = worldToCamera.rotation.cast<T>() * EigenVector3Map<T>(point3D) + worldToCamera.translation.cast<T>();
		CSimpleRadialCameraModel simpleRadialCameraModel;
		simpleRadialCameraModel.CameraToImage(cameraParams, point3DInCamera[0], point3DInCamera[1], point3DInCamera[2], &residuals[0], &residuals[1]);
		residuals[0] -= T(observedX);
		residuals[1] -= T(observedY);
		return true;
	}

private:
	const CRigid3D& worldToCamera;
	const double observedX;
	const double observedY;
};


inline void SetQuaternionManifold(ceres::Problem* problem, double* quatXYZW)
{
	problem->SetManifold(quatXYZW, new ceres::EigenQuaternionManifold());
}

inline void SetSubsetManifold(int size, const std::vector<int>& constantParams, ceres::Problem* problem, double* params)
{
	problem->SetManifold(params, new ceres::SubsetManifold(size, constantParams));
}

template <int size>
inline void SetSphereManifold(ceres::Problem* problem, double* params)
{
	problem->SetManifold(params, new ceres::SphereManifold<size>());
}





