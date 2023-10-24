#pragma once
#include <vector>
#include <flann/flann.hpp>

#include "Common.h"
#include "../Base/Options.h"
#include "../Scene/Point2D.h"
#include "../Scene/TwoViewGeometry.h"
#include "../../thirdparty/SIFT_GPU/SiftGPU.h"


class CSIFTCPUMatcher final
{
public:
	explicit CSIFTCPUMatcher(const CSIFTMatchingOptions& options);
	CSIFTMatches Match(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2);
	void MatchGuided(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, CTwoViewGeometry& twoViewGeometry);
	

private:
	CSIFTMatchingOptions options;
	using FlannIndexType = flann::Index<flann::L2<uint8_t>>;

	FlannIndexType BuildFlannIndex(const CSIFTDescriptors& descriptors);
};

class CSIFTGPUMatcher final
{
public:
	explicit CSIFTGPUMatcher(const CSIFTMatchingOptions& options);
	CSIFTMatches Match(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2);
	void MatchGuided(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, CTwoViewGeometry& twoViewGeometry);

	static size_t lastDescriptors1Index, lastDescriptors2Index;
	static bool isUploadDescriptors1, isUploadDescriptors2;
private:
	CSIFTMatchingOptions options;
	SiftMatchGPU siftMatchGPU;

};

