#include "Workflow.h"
#include "exif.h"
#include "../Base/Log.h"

using namespace std;

thread_local CSIFTCPUExtractor* CWorkflow::SIFTCPUExtractor = nullptr;
thread_local CSIFTGPUExtractor* CWorkflow::SIFTGPUExtractor = nullptr;
thread_local CSIFTCPUMatcher* CWorkflow::SIFTCPUMatcher = nullptr;
thread_local CSIFTGPUMatcher* CWorkflow::SIFTGPUMatcher = nullptr;

bool CWorkflow::ReadImage(const string& imagePath, CDatabase& database)
{
	auto start = chrono::high_resolution_clock::now();
	if (!IsFileExists(imagePath))
	{
		CLog::Log((boost::format("[Read image] The path %1% dose not exist!") % imagePath).str());
		return false;
	}

	double focalLength = 0;
	size_t width = 0, height = 0;
	GetExifData(imagePath, focalLength, width, height);
	CCamera camera;
	camera.SetCameraModelName("Simple Radial");
	if (focalLength > 0 && width != 0 && height != 0)
	{
		camera.SetFocalLength(focalLength);
		camera.SetFocalLengthPrior(true);
		camera.SetWidth(width);
		camera.SetHeight(height);
		camera.SetPrincipalPointX(width / 2.0);
		camera.SetPrincipalPointY(height / 2.0);
	}
	else if (width != 0 && height != 0)
	{
		camera.SetFocalLength(1.2 * max(width, height));
		camera.SetFocalLengthPrior(false);
		camera.SetWidth(width);
		camera.SetHeight(height);
		camera.SetPrincipalPointX(width / 2.0);
		camera.SetPrincipalPointY(height / 2.0);
	}
	else
	{
		const cv::Mat imageTemp = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
		if (imageTemp.empty() || imageTemp.cols <= 0 || imageTemp.rows <= 0)
		{
			CLog::Log((boost::format("[Read image] Error reading the image %1% !") % imagePath).str());
			return false;
		}
		camera.SetFocalLength(1.2 * max(imageTemp.cols, imageTemp.rows));
		camera.SetFocalLengthPrior(false);
		camera.SetWidth(imageTemp.cols);
		camera.SetHeight(imageTemp.rows);
		camera.SetPrincipalPointX(imageTemp.cols / 2.0);
		camera.SetPrincipalPointY(imageTemp.rows / 2.0);
	}
	Check(camera.VerifyParamsNum());

	size_t cameraID = numeric_limits<size_t>::max();
	if (!database.IsCameraExists(camera))
	{
		cameraID = database.AddCamera(camera);
	}
	else
	{
		cameraID = database.GetCameraID(camera);
	}
	CImage image;
	image.SetCameraID(cameraID);
	const string imageName = GetFileName(imagePath);
	image.SetImageName(imageName);
	const size_t imageID = database.AddImage(image);
	auto end = chrono::high_resolution_clock::now();
	size_t timeConsuming = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	CLog::Log((boost::format("[Read image %1%ms] Image %2% has been successfully read!") % to_string(timeConsuming) % imageName).str());
	return true;
}
bool CWorkflow::SIFTExtract(const string& imagePath, const CSIFTExtractionOptions& options, CDatabase& database)
{
	auto start = chrono::high_resolution_clock::now();
	const string imageName = GetFileName(imagePath);
	CKeypoints keypoints(0);
	CSIFTDescriptors descriptors;
	bool isSuccess = false;
	if (options.isUseGPU)
	{
		if (!SIFTGPUExtractor)
		{
			SIFTGPUExtractor = new CSIFTGPUExtractor(options);
		}
		isSuccess = SIFTGPUExtractor->Extract(imagePath, keypoints, descriptors);
	}
	else
	{
		if (!SIFTCPUExtractor)
		{
			SIFTCPUExtractor = new CSIFTCPUExtractor(options);
		}
		isSuccess = SIFTCPUExtractor->Extract(imagePath, keypoints, descriptors);
	}
	if (!isSuccess)
	{
		CLog::Log((boost::format("[Extraction] Error occurred while extracting SIFT features from image %1% !") % imageName).str());
		return false;
	}
	size_t numKeypoints = keypoints.size();
	Check(numKeypoints == descriptors.rows() && descriptors.cols() == 128);

	auto startWait = chrono::high_resolution_clock::now();
	while (!database.IsImageExists(imageName))
	{
		this_thread::sleep_for(chrono::milliseconds(50));
		if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startWait).count() >= 5)
		{
			CLog::Log((boost::format("[Extraction] Feature extraction timed out for Image %1% !") % imageName).str());
			return false;
		}
	}
	const size_t imageID = database.GetImageID(imageName);
	database.AddImageKeypoints(imageID, keypoints);
	database.AddImageDescriptors(imageID, descriptors);

	auto end = chrono::high_resolution_clock::now();
	size_t timeConsuming = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	CLog::Log((boost::format("[Extraction %1%ms] Image %2% had %3% feature points extracted.") % to_string(timeConsuming) % imageName % to_string(numKeypoints)).str());
	return true;
}
bool CWorkflow::SIFTMatch(CDatabase& database, size_t imageID1, size_t imageID2, const CSIFTMatchingOptions& options)
{
	if (!database.IsImageExists(imageID1) || !database.IsImageExists(imageID2))
	{
		return false;
	}
	if (database.GetImageDescriptorsNum(imageID1) == 0 || database.GetImageDescriptorsNum(imageID2) == 0)
	{
		return false;
	}
	const CSIFTDescriptors& descriptors1 = database.GetImageDescriptors(imageID1);
	const CSIFTDescriptors& descriptors2 = database.GetImageDescriptors(imageID2);
	Check(descriptors1.cols() == 128 && descriptors2.cols() == 128);

	CSIFTMatches matches;
	if (options.isUseGPU)
	{
		if (!SIFTGPUMatcher)
		{
			SIFTGPUMatcher = new CSIFTGPUMatcher(options);
		}

		CSIFTGPUMatcher::isUploadDescriptors1 = (CSIFTGPUMatcher::lastDescriptors1Index != imageID1);
		CSIFTGPUMatcher::isUploadDescriptors2 = (CSIFTGPUMatcher::lastDescriptors2Index != imageID2);
		matches = SIFTGPUMatcher->Match(descriptors1, descriptors2);
		CSIFTGPUMatcher::lastDescriptors1Index = imageID1;
		CSIFTGPUMatcher::lastDescriptors2Index = imageID2;
	}
	else
	{
		if (!SIFTCPUMatcher)
		{
			SIFTCPUMatcher = new CSIFTCPUMatcher(options);
		}
		matches = SIFTCPUMatcher->Match(descriptors1, descriptors2);
	}
	if (!matches.empty())
	{
		database.AddMatches(imageID1, imageID2, matches);
	}
	return true;
}
bool CWorkflow::EstimateTwoViewGeometry(CDatabase& database, size_t imageID1, size_t imageID2, const CTwoViewGeometryOptions& options)
{
	if (!database.IsImageExists(imageID1) || !database.IsImageExists(imageID2))
	{
		return false;
	}

	const CKeypoints keypoints1 = database.GetImageKeypoints(imageID1);
	const CKeypoints keypoints2 = database.GetImageKeypoints(imageID2);
	if (keypoints1.empty() || keypoints2.empty() || database.GetMatchesNum(imageID1, imageID2) == 0)
	{
		return true;
	}
	vector<Eigen::Vector2d> pointsVector1(keypoints1.size());
	vector<Eigen::Vector2d> pointsVector2(keypoints2.size());
	for (size_t i = 0; i < keypoints1.size(); i++)
	{
		pointsVector1[i].x() = keypoints1[i].pt.x;
		pointsVector1[i].y() = keypoints1[i].pt.y;
	}
	for (size_t i = 0; i < keypoints2.size(); i++)
	{
		pointsVector2[i].x() = keypoints2[i].pt.x;
		pointsVector2[i].y() = keypoints2[i].pt.y;
	}
	const CSIFTMatches matches = database.GetMatches(imageID1, imageID2);

	CTwoViewGeometry twoViewGeometry;
	twoViewGeometry.Estimate(database.GetImageCamera(imageID1), pointsVector1, database.GetImageCamera(imageID2), pointsVector2, matches, options);
	if (!twoViewGeometry.inlierMatches.empty())
	{
		database.AddTwoViewGeometry(imageID1, imageID2, twoViewGeometry);
	}
	return true;
}







