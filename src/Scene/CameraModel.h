#pragma once
#include <string>
#include <vector>
#include <limits>
#include <Eigen/Core>
#include <Eigen/LU>

#include "../Base/Base.h"

enum class CCameraModelType : size_t
{
	CSimplePinholeCameraModel = 0,			// 最基本的针孔相机模型, 不考虑镜头畸变. 参数: 焦距, 主点坐标
	CPinholeCameraModel = 1,				// 基本针孔模型的扩展, 包含畸变参数.
	CSimpleRadialCameraModel = 2,			// 考虑简单径向畸变的基本针孔模型. 参数: 焦距, 主点坐标, 径向畸变参数
	CSimpleRadialFisheyeCameraModel = 3,	// 针对鱼眼相机镜头的SimpleRadial模型
	CRadialCameraModel = 4,                 // 包含更复杂的径向畸变模型
	CRadialFisheyeCameraModel = 5,          // 针对鱼眼相机镜头的Radial模型
	COpenCVCameraModel = 6,                 // 以OpenCV的方式描述相机模型
	COpenCVFisheyeCameraModel = 7,          // 针对鱼眼相机镜头的OpenCV模型
	CFullOpenCVCameraModel = 8,             // 更全面和详细的OpenCV相机模型, 包含更多的参数和功能
	CFOVCameraModel = 9,                    // 用一个"视场"参数来描述径向畸变. 参数包括视场参数
	CThinPrismFisheyeCameraModel = 10,      // 结合了鱼眼镜头模型和薄棱镜畸变模型
	CCameraModelTypeNum = 11                // 共有11种相机模型
};

class CCameraModel
{
public:
	virtual std::string GetModelName() const noexcept;
	virtual size_t GetModelID() const noexcept;
	virtual std::string GetParamsType() const noexcept;
	virtual size_t GetParamsNum() const noexcept;
	virtual void GetFocalLengthIndex(size_t& focalLengthIndex_X, size_t& focalLengthIndex_Y) const noexcept;
	virtual void GetPrincipalPointIndex(size_t& principalPointIndex_X, size_t& principalPointIndex_Y) const noexcept;
	virtual std::vector<size_t> GetExtraParams() const noexcept;

	virtual bool IsBogusParams(const std::vector<double>& params, size_t width, size_t height, double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam) const;
	virtual bool IsBogusFocalLength(const std::vector<double>& params, size_t width, size_t height, double minFocalLengthRatio, double maxFocalLengthRatio) const;
	virtual bool IsBogusPrincipalPoint(const std::vector<double>& params, size_t width, size_t height) const;
	virtual bool IsBogusExtraParams(const std::vector<double>& params, double maxExtraParam) const;

	// 阈值除以焦距分量的平均值
	virtual double ImageToCameraThreshold(const std::vector<double>& params, double threshold) const;

	//初始化相机模型参数. 输入焦距, 影像宽, 影像高. 输出相机参数
	virtual std::vector<double> InitializeParams(double focalLength, size_t width, size_t height) const = 0;

	// 畸变引入: 从理想的相机坐标系坐标转到带畸变的影像坐标系坐标. w为归一化因子, 输出x和y为图像坐标系坐标
	virtual void CameraToImage(const std::vector<double>& params, float u, float v, float w, float& x, float& y) const = 0;

	// 畸变矫正: 从畸变的影像坐标系坐标转到理想的相机坐标系坐标.
	virtual void ImageToCamera(const std::vector<double>& params, float x, float y, float& u, float& v, float& w) const = 0;

	// 去畸变. 返回的u和v的新值在畸变之后会映射到u和v的原值
	virtual void IterativeUndistortion(const std::vector<double>& extraParams, float& u, float& v) const;

protected:
	std::string modelName;
	std::string paramsType;

	size_t modelID = std::numeric_limits<size_t>::max();
	size_t numParams = std::numeric_limits<size_t>::max();
	
	size_t focalLengthIndex_X = std::numeric_limits<size_t>::max(); //params的第focalLengthIndex_X处元素表示fx
	size_t focalLengthIndex_Y = std::numeric_limits<size_t>::max(); //params的第focalLengthIndex_Y处元素表示fy
	size_t principalPointIndex_X = std::numeric_limits<size_t>::max(); //params的第principalPointIndex_X处元素表示cx
	size_t principalPointIndex_Y = std::numeric_limits<size_t>::max(); //params的第principalPointIndex_Y处元素表示cy
	std::vector<size_t> extraParamsIndex;

	virtual const std::vector<double> GetExtraParams(const std::vector<double>& params) const = 0;
	virtual void Distortion(const std::vector<double>& extraParams, double u, double v, double& du, double& dv) const = 0;
};

// Simple Radial相机模型. 参数列表: f, cx, cy, k
class CSimpleRadialCameraModel final :public CCameraModel
{
public:
	CSimpleRadialCameraModel() noexcept;
	std::vector<double> InitializeParams(double focalLength, size_t width, size_t height) const noexcept override;
	void CameraToImage(const std::vector<double>& params, float u, float v, float w, float& x, float& y) const override;

	template <typename T>
	static void CameraToImage(const T* params, T u, T v, T w, T* x, T* y);
	void ImageToCamera(const std::vector<double>& params, float x, float y, float& u, float& v, float& w) const override;

private:
	const std::vector<double> GetExtraParams(const std::vector<double>& params) const override;
	void Distortion(const std::vector<double>& extraParams, double u, double v, double& du, double& dv) const override;

	template <typename T>
	static void Distortion(const T* extraParams, T u, T v, T* du, T* dv);
};












