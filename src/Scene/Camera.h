#pragma once
#include <vector>
#include <string>

#include "CameraModel.h"

// ���. �����ڲ�, �����ж��ͼ�������(����: ͬһ���������ʹ����ȫ��ͬ�ľ�ͷ���ڲ������˶���Ӱ��)
// ������������ģ���ඨ����ض�����ģ��
class CCamera final
{
public:
	CCamera() noexcept;
	~CCamera() = default;

	// ͨ���������ģ��ID�ͽ����ʼ�����������������������Ϊͼ�������
	CCamera(size_t modelID, double focalLength, size_t width, size_t height);

	// ͨ���������ģ�����ƺͽ����ʼ�����������������������Ϊͼ�������
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
	
	// ��ȡ�ڲξ���(�������������, �������������)
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

	// ��ȡ��������
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

	// ��һ���������ŷָ��Ĳ������ַ�����������Ĳ���
	bool SetParamsFromString(const std::string& paramsString);
	inline const double* GetParamsData() const noexcept
	{
		return params.data();
	}
	inline double* GetParamsData() noexcept
	{
		return params.data();
	}

	// ����������Զ��ŷָ����ַ�����ʽ����
	std::string GetParamsString() const;

	// �������Ĳ����Ƿ���Ч��������������ά���Ƿ���ָ�������ģ��ƥ��
	inline bool VerifyParamsNum() const noexcept
	{
		if (!cameraModel)
		{
			return params.empty();
		}
		return (params.size() == cameraModel->GetParamsNum());
	}

	// �������Ƿ���в�����Ĳ���
	inline bool IsBogusParams(double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam) const
	{
		if (!cameraModel)
		{
			return false;
		}
		return cameraModel->IsBogusParams(params, width, height, minFocalLengthRatio, maxFocalLengthRatio, maxExtraParam);
	}

	// �������Ƿ��޻���
	bool IsUndistorted() const;

	// ����У��: �ӻ����Ӱ������ϵ����ת��������������ϵ����
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

	// ��Ӱ��ƽ���е���ֵת�����������ϵ�µ���ֵ
	inline double ImageToCameraThreshold(double threshold) const
	{
		CHECK(cameraModel != nullptr);
		return cameraModel->ImageToCameraThreshold(params, threshold);
	}
	
	// ��������: ��������������ϵ����ת���������Ӱ������ϵ����
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

	// �������ά��, �������������
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
	size_t ID = std::numeric_limits<size_t>::max();  // ���û��ָ�����ID, ��ᱻ��Ϊstd::numeric_limits<size_t>::max()
	std::shared_ptr<CCameraModel> cameraModel = nullptr;  // ���ģ��(����ָ��ͬ������)
	size_t width, height;                            // Ӱ��ߴ�
	std::vector<double> params;                      // ���ࡢ����Ͷ������. ���û��ָ�����ģ��, ��ô��vectorΪ��
	bool isFocalLengthPrior = false;                 // �Ƿ��пɿ������齹����Ϣ(����: �ֶ���ȡ���ߴ�Ӱ���EXIF��Ϣ����ȡ)
};










