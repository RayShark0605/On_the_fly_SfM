#include "CameraModel.h"

using namespace std;

string CCameraModel::GetModelName() const noexcept
{
	return modelName;
}
size_t CCameraModel::GetModelID() const noexcept
{
	return modelID;
}
string CCameraModel::GetParamsType() const noexcept
{
	return paramsType;
}
size_t CCameraModel::GetParamsNum() const noexcept
{
	return numParams;
}
void CCameraModel::GetFocalLengthIndex(size_t& focalLengthIndex_X, size_t& focalLengthIndex_Y) const noexcept
{
	focalLengthIndex_X = this->focalLengthIndex_X;
	focalLengthIndex_Y = this->focalLengthIndex_Y;
}
void CCameraModel::GetPrincipalPointIndex(size_t& principalPointIndex_X, size_t& principalPointIndex_Y) const noexcept
{
	principalPointIndex_X = this->principalPointIndex_X;
	principalPointIndex_Y = this->principalPointIndex_Y;
}
vector<size_t> CCameraModel::GetExtraParams() const noexcept
{
	return extraParamsIndex;
}
bool CCameraModel::IsBogusParams(const vector<double>& params, size_t width, size_t height, double minFocalLengthRatio, double maxFocalLengthRatio, double maxExtraParam) const
{
	return IsBogusPrincipalPoint(params, width, height) ||
		IsBogusFocalLength(params, width, height, minFocalLengthRatio, maxExtraParam) ||
		IsBogusExtraParams(params, maxExtraParam);
}
bool CCameraModel::IsBogusFocalLength(const vector<double>& params, size_t width, size_t height, double minFocalLengthRatio, double maxFocalLengthRatio) const
{
	Check(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());
	const double invMaxSize = 1.0 / max(width, height);
	const double invFocalLengthX = params[focalLengthIndex_X] * invMaxSize;
	const double invFocalLengthY = params[focalLengthIndex_Y] * invMaxSize;
	return (invFocalLengthX < minFocalLengthRatio || invFocalLengthX > maxFocalLengthRatio ||
		invFocalLengthY < minFocalLengthRatio || invFocalLengthY> maxFocalLengthRatio);
}
bool CCameraModel::IsBogusPrincipalPoint(const vector<double>& params, size_t width, size_t height) const
{
	Check(principalPointIndex_X < params.size() && principalPointIndex_Y < params.size());
	const double cx = params[principalPointIndex_X];
	const double cy = params[principalPointIndex_Y];
	return (cx < 0 || cx > width || cy < 0 || cy > height);
}
bool CCameraModel::IsBogusExtraParams(const vector<double>& params, double maxExtraParam) const
{
	for (const size_t index : extraParamsIndex)
	{
		Check(index < params.size());
		if (abs(params[index]) > maxExtraParam)
		{
			return true;
		}
	}
	return false;
}
double CCameraModel::ImageToCameraThreshold(const vector<double>& params, double threshold) const
{
	Check(focalLengthIndex_X < params.size() && focalLengthIndex_Y < params.size());
	double meanFocalLength = params[focalLengthIndex_X] / 2.0 + params[focalLengthIndex_Y] / 2.0;
	return threshold / meanFocalLength;
}
void CCameraModel::IterativeUndistortion(const std::vector<double>& extraParams, float& u, float& v) const
{
	const size_t numIterations = 100; // Newton迭代方法的最大迭代次数
	const double maxStepNorm = 1e-10; // 当迭代步长的范数(每次迭代的变化量)小于这个值时，迭代应该停止
	const double relStepSize = 1e-6;  // 计算中心差分(一种数值微分方法)的相对步长
	Eigen::Matrix2d J; // 雅可比矩阵
	const Eigen::Vector2d x0(u, v);
	Eigen::Vector2d x(u, v);
	Eigen::Vector2d dx;
	Eigen::Vector2d dx_0b;
	Eigen::Vector2d dx_0f;
	Eigen::Vector2d dx_1b;
	Eigen::Vector2d dx_1f;
	for (size_t i = 0; i < numIterations; i++)
	{
		const double step0 = std::max(std::numeric_limits<double>::epsilon(), std::abs(relStepSize * x(0)));
		const double step1 = std::max(std::numeric_limits<double>::epsilon(), std::abs(relStepSize * x(1)));
		Distortion(extraParams, x(0), x(1), dx(0), dx(1));
		Distortion(extraParams, x(0) - step0, x(1), dx_0b(0), dx_0b(1));
		Distortion(extraParams, x(0) + step0, x(1), dx_0f(0), dx_0f(1));
		Distortion(extraParams, x(0), x(1) - step1, dx_1b(0), dx_1b(1));
		Distortion(extraParams, x(0), x(1) + step1, dx_1f(0), dx_1f(1));
		J(0, 0) = 1 + (dx_0f(0) - dx_0b(0)) / (2 * step0);
		J(0, 1) = (dx_1f(0) - dx_1b(0)) / (2 * step1);
		J(1, 0) = (dx_0f(1) - dx_0b(1)) / (2 * step0);
		J(1, 1) = 1 + (dx_1f(1) - dx_1b(1)) / (2 * step1);
		const Eigen::Vector2d stepX = J.partialPivLu().solve(x + dx - x0);
		x -= stepX;
		if (stepX.squaredNorm() < maxStepNorm)
		{
			break;
		}
	}
	u = x(0);
	v = x(1);
}


CSimpleRadialCameraModel::CSimpleRadialCameraModel() noexcept
{
	modelName = "Simple Radial";
	paramsType = "f, cx, cy, k";
	modelID = static_cast<size_t>(CCameraModelType::CSimpleRadialCameraModel);
	numParams = 4;
	focalLengthIndex_X = focalLengthIndex_Y = 0;
	principalPointIndex_X = 1;
	principalPointIndex_Y = 2;
	extraParamsIndex = { 3 };
}
vector<double> CSimpleRadialCameraModel::InitializeParams(double focalLength, size_t width, size_t height) const noexcept
{
	return { focalLength, width / 2.0, height / 2.0, 0 };
}
void CSimpleRadialCameraModel::CameraToImage(const vector<double>& params, float u, float v, float w, float& x, float& y) const
{
	Check(focalLengthIndex_X == focalLengthIndex_Y);
	Check(focalLengthIndex_X < params.size());
	Check(principalPointIndex_X < params.size());
	Check(principalPointIndex_Y < params.size());

	// 提取相机参数
	const double f = params[focalLengthIndex_X];
	const double c1 = params[principalPointIndex_X];
	const double c2 = params[principalPointIndex_Y];

	// 归一化
	u /= w;
	v /= w;

	// 畸变引入
	double du, dv;
	Distortion(GetExtraParams(params), u, v, du, dv);
	x = u + du;
	y = v + dv;

	// 转换到图像坐标系
	x = f * x + c1;
	y = f * y + c2;
}
template <typename T>
void CSimpleRadialCameraModel::CameraToImage(const T* params, T u, T v, T w, T* x, T* y)
{
	const double f = params[0];
	const double c1 = params[1];
	const double c2 = params[2];
	u /= w;
	v /= w;
	// 畸变引入
	double du, dv;
	Distortion(&params[3], u, v, &du, &dv);
	*x = u + du;
	*y = v + dv;

	// 转换到图像坐标系
	*x = f * (*x) + c1;
	*y = f * (*y) + c2;

}
void CSimpleRadialCameraModel::ImageToCamera(const vector<double>& params, float x, float y, float& u, float& v, float& w) const
{
	Check(focalLengthIndex_X == focalLengthIndex_Y);
	Check(focalLengthIndex_X < params.size());
	Check(principalPointIndex_X < params.size());
	Check(principalPointIndex_Y < params.size());

	// 提取相机参数
	const float f = params[focalLengthIndex_X];
	const float c1 = params[principalPointIndex_X];
	const float c2 = params[principalPointIndex_Y];

	// 将图像坐标转到归一化平面, w为1表示是一个同质坐标
	u = (x - c1) / f;
	v = (y - c2) / f;
	w = 1;
	
	// 迭代去畸变
	IterativeUndistortion(GetExtraParams(params), u, v);
}
const vector<double> CSimpleRadialCameraModel::GetExtraParams(const vector<double>& params) const
{
	Check(params.size() > 3);
	return { params[3] };
}
void CSimpleRadialCameraModel::Distortion(const vector<double>& extraParams, double u, double v, double& du, double& dv) const
{
	Check(extraParams.size() == 1);

	const double k = extraParams[0];
	const double u2 = u * u;
	const double v2 = v * v;
	const double r2 = u2 + v2;
	const double radial = k * r2;
	du = u * radial;
	dv = v * radial;
}
template <typename T>
static void CSimpleRadialCameraModel::Distortion(const T* extraParams, T u, T v, T* du, T* dv)
{
	const T k = extraParams[0];

	const T u2 = u * u;
	const T v2 = v * v;
	const T r2 = u2 + v2;
	const T radial = k * r2;
	*du = u * radial;
	*dv = v * radial;
}