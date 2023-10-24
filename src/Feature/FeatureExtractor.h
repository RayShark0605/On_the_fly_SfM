#pragma once
#include <string>
#include <vector>

#include "../../src/Base/Options.h"
#include "../../src/Scene/Point2D.h"
#include "../../thirdparty/VLFeat/sift.h"
#include "../../thirdparty/SIFT_GPU/SiftGPU.h"
#include "Common.h"

class CSIFTExtractor
{
public:
	virtual bool Extract(const std::string& imagePath, CKeypoints& keypoints, CSIFTDescriptors& descriptors) = 0;
	~CSIFTExtractor() = default;
};


class CSIFTCPUExtractor final :public CSIFTExtractor
{
public:
	using VlSiftType = std::unique_ptr<VlSiftFilt, void (*)(VlSiftFilt*)>;
	explicit CSIFTCPUExtractor(const CSIFTExtractionOptions& options);
	bool Extract(const std::string& imagePath, CKeypoints& keypoints, CSIFTDescriptors& descriptors) override;

private:
	CSIFTExtractionOptions options;
	VlSiftType sift;

	bool Extract_VLFeat(const std::vector<uchar>& imageData, size_t width, size_t height, CKeypoints& keypoints, CSIFTDescriptors& descriptors);
	bool Extract_OpenCV(const cv::Mat& image, CKeypoints& keypoints, CSIFTDescriptors& descriptors);
};

class CSIFTGPUExtractor final :public CSIFTExtractor
{
public:
	explicit CSIFTGPUExtractor(const CSIFTExtractionOptions& options);
	bool Extract(const std::string& imagePath, CKeypoints& keypoints, CSIFTDescriptors& descriptors) override;

private:
	CSIFTExtractionOptions options;
	SiftGPU siftGPU;

	bool Extract_SIFTGPU(const std::vector<uchar>& imageData, size_t width, size_t height, CKeypoints& keypoints, CSIFTDescriptors& descriptors);
};















