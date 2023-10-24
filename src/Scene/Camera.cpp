#include "Camera.h"

using namespace std;


CCamera::CCamera() noexcept
{
	ID = numeric_limits<size_t>::max();
	cameraModel = nullptr;
	width = height = 0;
	params = {};
	isFocalLengthPrior = false;
}
CCamera::CCamera(size_t modelID, double focalLength, size_t width, size_t height)
{
	SetCameraModelID(modelID);
	params = cameraModel->InitializeParams(focalLength, width, height);
	this->width = width;
	this->height = height;
}
CCamera::CCamera(const string& modelName, double focalLength, size_t width, size_t height)
{
	SetCameraModelName(modelName);
	params = cameraModel->InitializeParams(focalLength, width, height);
	this->width = width;
	this->height = height;
}
string CCamera::GetParamsString() const
{
	return VectorToString(params);
}
bool CCamera::SetParamsFromString(const string& paramsString)
{
	if (!cameraModel)
	{
		return false;
	}
	const vector<double> newParams = StringToVector<double>(paramsString);
	if (newParams.size() != cameraModel->GetParamsNum())
	{
		return false;
	}
	params = newParams;
	return true;
}
bool CCamera::IsUndistorted() const
{
	vector<size_t> extraParamsIndex = GetExtraParamsIndex();
	for (const size_t index : extraParamsIndex)
	{
		CHECK(index < params.size());
		if (abs(params[index] > 1e-8))
		{
			return false;
		}
	}
	return true;
}
