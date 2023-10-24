#pragma once
#include <string>
#include "../Scene/Database.h"
#include "../Base/Options.h"
#include "../Base/Base.h"
#include "../Feature/FeatureExtractor.h"
#include "../Feature/FeatureMatcher.h"
#include "ThreadPool.h"


class CWorkflow final
{
public:
	static bool ReadImage(const std::string& imagePath, CDatabase& database);
	static bool SIFTExtract(const std::string& imagePath, const CSIFTExtractionOptions& options, CDatabase& database);
	static bool SIFTMatch(CDatabase& database, size_t imageID1, size_t imageID2, const CSIFTMatchingOptions& options);
	static bool EstimateTwoViewGeometry(CDatabase& database, size_t imageID1, size_t imageID2, const CTwoViewGeometryOptions& options);

private:
	static thread_local CSIFTCPUExtractor* SIFTCPUExtractor;
	static thread_local CSIFTGPUExtractor* SIFTGPUExtractor;
	static thread_local CSIFTCPUMatcher* SIFTCPUMatcher;
	static thread_local CSIFTGPUMatcher* SIFTGPUMatcher;


};







