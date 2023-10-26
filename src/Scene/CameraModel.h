#pragma once
#include <string>
#include <vector>
#include <limits>
#include <Eigen/Core>
#include <Eigen/LU>

#include "../Base/Base.h"

enum class CCameraModelType : size_t
{
	CSimplePinholeCameraModel = 0,			// �������������ģ��, �����Ǿ�ͷ����. ����: ����, ��������
	CPinholeCameraModel = 1,				// �������ģ�͵���չ, �����������.
	CSimpleRadialCameraModel = 2,			// ���Ǽ򵥾������Ļ������ģ��. ����: ����, ��������, ����������
	CSimpleRadialFisheyeCameraModel = 3,	// ������������ͷ��SimpleRadialģ��
	CRadialCameraModel = 4,                 // ���������ӵľ������ģ��
	CRadialFisheyeCameraModel = 5,          // ������������ͷ��Radialģ��
	COpenCVCameraModel = 6,                 // ��OpenCV�ķ�ʽ�������ģ��
	COpenCVFisheyeCameraModel = 7,          // ������������ͷ��OpenCVģ��
	CFullOpenCVCameraModel = 8,             // ��ȫ�����ϸ��OpenCV���ģ��, ��������Ĳ����͹���
	CFOVCameraModel = 9,                    // ��һ��"�ӳ�"�����������������. ���������ӳ�����
	CThinPrismFisheyeCameraModel = 10,      // ��������۾�ͷģ�ͺͱ��⾵����ģ��
	CCameraModelTypeNum = 11                // ����11�����ģ��
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

	// ��ֵ���Խ��������ƽ��ֵ
	virtual double ImageToCameraThreshold(const std::vector<double>& params, double threshold) const;

	//��ʼ�����ģ�Ͳ���. ���뽹��, Ӱ���, Ӱ���. ����������
	virtual std::vector<double> InitializeParams(double focalLength, size_t width, size_t height) const = 0;

	// ��������: ��������������ϵ����ת���������Ӱ������ϵ����. wΪ��һ������, ���x��yΪͼ������ϵ����
	virtual void CameraToImage(const std::vector<double>& params, float u, float v, float w, float& x, float& y) const = 0;

	// �������: �ӻ����Ӱ������ϵ����ת��������������ϵ����.
	virtual void ImageToCamera(const std::vector<double>& params, float x, float y, float& u, float& v, float& w) const = 0;

	// ȥ����. ���ص�u��v����ֵ�ڻ���֮���ӳ�䵽u��v��ԭֵ
	virtual void IterativeUndistortion(const std::vector<double>& extraParams, float& u, float& v) const;

protected:
	std::string modelName;
	std::string paramsType;

	size_t modelID = std::numeric_limits<size_t>::max();
	size_t numParams = std::numeric_limits<size_t>::max();
	
	size_t focalLengthIndex_X = std::numeric_limits<size_t>::max(); //params�ĵ�focalLengthIndex_X��Ԫ�ر�ʾfx
	size_t focalLengthIndex_Y = std::numeric_limits<size_t>::max(); //params�ĵ�focalLengthIndex_Y��Ԫ�ر�ʾfy
	size_t principalPointIndex_X = std::numeric_limits<size_t>::max(); //params�ĵ�principalPointIndex_X��Ԫ�ر�ʾcx
	size_t principalPointIndex_Y = std::numeric_limits<size_t>::max(); //params�ĵ�principalPointIndex_Y��Ԫ�ر�ʾcy
	std::vector<size_t> extraParamsIndex;

	virtual const std::vector<double> GetExtraParams(const std::vector<double>& params) const = 0;
	virtual void Distortion(const std::vector<double>& extraParams, double u, double v, double& du, double& dv) const = 0;
};

// Simple Radial���ģ��. �����б�: f, cx, cy, k
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












