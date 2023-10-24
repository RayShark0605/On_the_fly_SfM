#pragma once
#include <vector>
#include <string>

#include "CameraModel.h"

// 相机. 包含内参, 可能有多个图像共享相机(例如: 同一个物理相机使用完全相同的镜头和内参拍摄了多张影像)
// 该类具有由相机模型类定义的特定畸变模型
class CCamera final
{
public:
	CCamera() noexcept;
	~CCamera() = default;

	// 通过给定相机模型ID和焦距初始化相机参数，并将主点设置为图像的中心
	CCamera(size_t modelID, double focalLength, size_t width, size_t height);

	// 通过给定相机模型名称和焦距初始化相机参数，并将主点设置为图像的中心
	CCamera(const std::string& modelName, double focalLength, size_t width, size_t height);

	bool operator==(const CCamera& other) const
	{
		bool isCameraModelEqual = (cameraModel == nullptr || other.cameraModel == nullptr) ? true : (cameraModel->GetModelID() == other.cameraModel->GetModelID());
		return (isCameraModelEqual && width == other.width && height == other.height && isFocalLengthPrior == other.isFocalLengthPrior && params == other.params);
	}

	inline size_t GetCameraID() const noexcept
	{
		return ID;
	}
	inline void SetCameraID(size_t cameraID) noexcept
	{
		ID = cameraID;
	}
	inline size_t GetCameraModelID() const noexcept
	{
		if (!cameraModel)
		{
			return std::numeric_limits<size_t>::max();
		}
		return cameraModel->GetModelID();
	}
	inline std::string GetCameraModelName() const noexcept
	{
		if (!cameraModel)
		{
			return "Unknown";
		}
		return cameraModel->GetModelName();
	}
	inline void SetCameraModelID(size_t modelID)
	{
		switch (modelID)
		{
		case static_cast<size_t>(CCameraModelType::CSimpleRadialCameraModel):
			cameraModel = std::make_shared<CSimpleRadialCameraModel>();
			params.resize(cameraModel->GetParamsNum());
			break;
		default:
			CHECK(false, "Unknown camera model ID!");
		}
	}
	inline void SetCameraModelName(const std::string& modelName)
	{
		if (modelName == "Simple Radial")
		{
			cameraModel = std::make_shared<CSimpleRadialCameraModel>();
			params.resize(cameraModel->GetParamsNum());
		}
		else
		{
			CHECK(false, "Unknown camera model name!");
		}
	}

	inline size_t GetWidth() const noexcept
	{
		return width;
	}
	inline size_t GetHeight() const noexcept
	{
		return height;
	}
	inline void SetWidth(size_t width) noexcept
	{
		this->width = width;
	}
	inline void SetHeight(size_t height) noexcept
	{
		this->height = height;
	}

	inline void GetFocalLengthIndex(size_t& focalLengthIndex_X, size_t& focalLengthIndex_Y) const noexcept
	{
		if (!cameraModel)
		{
			focalLengthIndex_X = focalLengthIndex_Y = std::numeric_limits<size_t>::max();
			return;
		}
		cameraModel->GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
	}
	inline bool IsFocalLengthPrior() const noexcept
	{
		return isFocalLengthPrior;
	}
	inline void SetFocalLengthPrior(bool isFocalLengthPrior) noexcept
	{
		this->isFocalLengthPrior = isFocalLengthPrior;
	}
	inline double GetMeanFocalLength() const
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());

		double focalLength = params[focalLengthIndex_X] + params[focalLengthIndex_Y];
		return focalLength / 2.0;
	}
	inline double GetFocalLength() const
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);

		CHECK(focalLengthIndex_X == focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size());
		return params[focalLengthIndex_X];
	}
	inline double GetFocalLengthX() const
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());
		return params[focalLengthIndex_X];
	}
	inline double GetFocalLengthY() const
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());
		return params[focalLengthIndex_Y];
	}
	inline void SetFocalLength(double focalLength)
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());

		params[focalLengthIndex_X] = focalLength;
		params[focalLengthIndex_Y] = focalLength;
	}
	inline void SetFocalLengthX(double focalLengthX)
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X != focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());

		params[focalLengthIndex_X] = focalLengthX;
	}
	inline void SetFocalLengthY(double focalLengthY)
	{
		size_t focalLengthIndex_X, focalLengthIndex_Y;
		GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X != focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());

		params[focalLengthIndex_Y] = focalLengthY;
	}

	inline void GetPrincipalPointIndex(size_t& principalPointIndex_X, size_t& principalPointIndex_Y) const noexcept
	{
		if (!cameraModel)
		{
			principalPointIndex_X = principalPointIndex_Y = std::numeric_limits<size_t>::max();
			return;
		}
		cameraModel->GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
	}
	inline double GetPrincipalPointX() const
	{
		size_t principalPointIndex_X, principalPointIndex_Y;
		GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
		CHECK(principalPointIndex_X < params.size() && principalPointIndex_Y < params.size());

		return params[principalPointIndex_X];
	}
	inline double GetPrincipalPointY() const
	{
		size_t principalPointIndex_X, principalPointIndex_Y;
		GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
		CHECK(principalPointIndex_X < params.size() && principalPointIndex_Y < params.size());

		return params[principalPointIndex_Y];
	}
	inline void SetPrincipalPointX(double principalPointX)
	{
		size_t principalPointIndex_X, principalPointIndex_Y;
		GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
		CHECK(principalPointIndex_X != principalPointIndex_Y);
		CHECK(principalPointIndex_X < params.size() && principalPointIndex_Y < params.size());

		params[principalPointIndex_X] = principalPointX;
	}
	inline void SetPrincipalPointY(double principalPointY)
	{
		size_t principalPointIndex_X, principalPointIndex_Y;
		GetPrincipalPointIndex(principalPointIndex_X, principalPointIndex_Y);
		CHECK(principalPointIndex_X != principalPointIndex_Y);
		CHECK(principalPointIndex_X < params.size() && principalPointIndex_Y < params.size());

		params[principalPointIndex_Y] = principalPointY;
	}
	
	// 获取内参矩阵(包括焦距和主点, 不包括畸变参数)
	inline Eigen::Matrix3d GetCalibrationMatrix() const
	{
		CHECK(cameraModel != nullptr);

		Eigen::Matrix3d K = Eigen::Matrix3d::Identity();

		size_t focalLengthIndex_X, focalLengthIndex_Y;
		cameraModel->GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		CHECK(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());
		K(0, 0) = params[focalLengthIndex_X];
		K(1, 1) = params[focalLengthIndex_Y];
		K(0, 2) = GetPrincipalPointX();
		K(1, 2) = GetPrincipalPointY();
		return K;
	}

	// 获取参数类型
	inline std::string GetParamsType() const
	{
		CHECK(cameraModel != nullptr);
		return cameraModel->GetParamsType();
	}
	inline size_t GetParamsNum() const
	{
		CHECK(cameraModel != nullptr);
		return cameraModel->GetParamsNum();
	}
	inline const std::vector<double>& GetParams() const noexcept
	{
		return params;
	}
	inline std::vector<double>& GetParams() noexcept
	{
		return params;
	}
	inline double GetParams(size_t index) const
	{
		CHECK(index < params.size());
		return params[index];
	}
	inline double& GetParams(size_t index)
	{
		CHECK(index < params.size());
		return params[index];
	}
	inline void SetParams(const std::vector<double>& newParams) noexcept
	{
		params = newParams;
	}
	inline std::vector<size_t> GetExtraParamsIndex() const noexcept
	{
		if (!cameraModel)
		{
			return {};
		}
		return cameraModel->GetExtraParams();
	}

	// 从一个包含逗号分隔的参数的字符串设置相机的参数
	bool SetParamsFromString(const std::string& paramsString);
	inline const double* GetParamsData() const noexcept
	{
		return params.data();
	}
	inline double* GetParamsData() noexcept
	{
		return params.data();
	}

	// 将相机参数以逗号分隔的字符串形式返回
	std::string GetParamsString() const;

	// 检查相机的参数是否有效，即参数向量的维度是否与指定的相机模型匹配
	inline bool VerifyParamsNum() const noexcept
	{
		if (!cameraModel)
		{
			return params.empty();
		}
		return (params.size() == cameraModel->GetParamsNum());
	}

	// 检查相机是否具有不合理的参数
	inline bool IsBogusParams(double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam) const
	{
		if (!cameraModel)
		{
			return false;
		}
		return cameraModel->IsBogusParams(params, width, height, minFocalLengthRatio, maxFocalLengthRatio, maxExtraParam);
	}

	// 检查相机是否无畸变
	bool IsUndistorted() const;

	// 畸变校正: 从畸变的影像坐标系坐标转到理想的相机坐标系坐标
	inline void ImageToCamera(float imagePointX, float imagePointY, float& cameraPointX, float& cameraPointY) const
	{
		CHECK(cameraModel != nullptr);
		float u, v, w;
		cameraModel->ImageToCamera(params, imagePointX, imagePointY, u, v, w);
		CHECK(abs(w) > 1e-6);
		cameraPointX = u / w;
		cameraPointY = v / w;
	}
	inline Eigen::Vector2d ImageToCamera(const Eigen::Vector2d& imagePoint) const
	{
		CHECK(cameraModel != nullptr);
		float u, v, w;
		cameraModel->ImageToCamera(params, imagePoint(0), imagePoint(1), u, v, w);
		CHECK(abs(w) > 1e-6);
		return Eigen::Vector2d(u / w, v / w);
	}

	// 将影像平面中的阈值转换到相机坐标系下的阈值
	inline double ImageToCameraThreshold(double threshold) const
	{
		CHECK(cameraModel != nullptr);
		return cameraModel->ImageToCameraThreshold(params, threshold);
	}
	
	// 畸变引入: 从理想的相机坐标系坐标转到带畸变的影像坐标系坐标
	inline void CameraToImage(float cameraPointX, float cameraPointY, float& imagePointX, float& imagePointY) const
	{
		CHECK(cameraModel != nullptr);
		cameraModel->CameraToImage(params, cameraPointX, cameraPointY, 1, imagePointX, imagePointY);
	}
	inline Eigen::Vector2d CameraToImage(const Eigen::Vector2d& cameraPoint) const
	{
		CHECK(cameraModel != nullptr);
		float imagePointX, imagePointY;
		cameraModel->CameraToImage(params, cameraPoint(0), cameraPoint(1), 1, imagePointX, imagePointY);
		return Eigen::Vector2d(imagePointX, imagePointY);
	}

	// 缩放相机维度, 调整焦距和主点
	inline void Rescale(double scale)
	{
		CHECK(scale > 0);
		CHECK(cameraModel != nullptr);

		const double scaleX = std::round(scale * width) / width;
		const double scaleY = std::round(scale * height) / height;
		width = static_cast<size_t>(std::round(scale * width));
		height = static_cast<size_t>(std::round(scale * height));

		SetPrincipalPointX(scaleX * GetPrincipalPointX());
		SetPrincipalPointY(scaleY * GetPrincipalPointY());

		size_t focalLengthIndex_X, focalLengthIndex_Y;
		cameraModel->GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		if (focalLengthIndex_X == focalLengthIndex_Y)
		{
			SetFocalLength((scaleX + scaleY) / 2.0 * GetFocalLength());
		}
		else
		{
			SetFocalLengthX(scaleX * GetFocalLengthX());
			SetFocalLengthY(scaleY * GetFocalLengthY());
		}
	}
	inline void Rescale(size_t newWidth, size_t newHeight)
	{
		CHECK(cameraModel != nullptr);

		const double scaleX = newWidth * 1.0 / width;
		const double scaleY = newHeight * 1.0 / height;

		width = newWidth;
		height = newHeight;

		SetPrincipalPointX(scaleX * GetPrincipalPointX());
		SetPrincipalPointY(scaleY * GetPrincipalPointY());

		size_t focalLengthIndex_X, focalLengthIndex_Y;
		cameraModel->GetFocalLengthIndex(focalLengthIndex_X, focalLengthIndex_Y);
		if (focalLengthIndex_X == focalLengthIndex_Y)
		{
			SetFocalLength((scaleX + scaleY) / 2.0 * GetFocalLength());
		}
		else
		{
			SetFocalLengthX(scaleX * GetFocalLengthX());
			SetFocalLengthY(scaleY * GetFocalLengthY());
		}
	}

private:
	size_t ID = std::numeric_limits<size_t>::max();  // 如果没有指定相机ID, 则会被设为std::numeric_limits<size_t>::max()
	std::shared_ptr<CCameraModel> cameraModel = nullptr;  // 相机模型(可能指向不同的子类)
	size_t width, height;                            // 影像尺寸
	std::vector<double> params;                      // 焦距、主点和额外参数. 如果没有指定相机模型, 那么此vector为空
	bool isFocalLengthPrior = false;                 // 是否有可靠的先验焦距信息(例如: 手动提取或者从影像的EXIF信息中提取)
};










