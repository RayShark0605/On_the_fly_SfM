#include "FeatureMatcher.h"

using namespace std;

CSIFTCPUMatcher::CSIFTCPUMatcher(const CSIFTMatchingOptions& options)
{
	options.Check();
	this->options = options;
}
CSIFTMatches CSIFTCPUMatcher::Match(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2)
{
	CHECK(descriptors1.cols() == 128 && descriptors2.cols() == 128);
	CHECK(descriptors1.rows() != 0 && descriptors2.rows() != 0);

	FlannIndexType flannIndex1 = BuildFlannIndex(descriptors1);
	FlannIndexType flannIndex2 = BuildFlannIndex(descriptors2);

	CSIFTMatches matches;
	if (!options.isUseFLANN)
	{
		const Eigen::MatrixXi distances = ComputeSIFTDistanceMatrix(descriptors1, descriptors2);
		FindBestMatchesBruteForce(distances, options.maxRatio, options.maxDistance, options.doCrossCheck, matches);
		return matches;
	}
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices_1to2;
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> distances_1to2;

	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices_2to1;
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> distances_2to1;
	FindNearestNeighborsFlann(descriptors1, descriptors2, flannIndex2, indices_1to2, distances_1to2);
	if (options.doCrossCheck)
	{
		FindNearestNeighborsFlann(descriptors2, descriptors1, flannIndex1, indices_2to1, distances_2to1);
	}
	matches = FindBestMatchesFlann(indices_1to2, distances_1to2, indices_2to1, distances_2to1, options.maxRatio, options.maxDistance, options.doCrossCheck);
	return matches;
}
void CSIFTCPUMatcher::MatchGuided(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, CTwoViewGeometry& twoViewGeometry)
{
	CHECK(descriptors1.cols() == 128 && descriptors2.cols() == 128);
	CHECK(descriptors1.rows() != 0 && descriptors2.rows() != 0);
	CHECK(keypoints1.size() == descriptors1.rows() && keypoints2.size() == descriptors2.rows());
	if (twoViewGeometry.type < 2 || twoViewGeometry.type > 6)
	{
		return;
	}

	float max_residual = options.maxError * options.maxError;
	const Eigen::Matrix3f F = twoViewGeometry.F.cast<float>();
	const Eigen::Matrix3f H = twoViewGeometry.H.cast<float>();
	size_t inlierMatchesNum = twoViewGeometry.inlierMatches.size();
	CSIFTMatches newInlierMatches;
	newInlierMatches.reserve(inlierMatchesNum);
	for (size_t i = 0; i < inlierMatchesNum; i++)
	{
		const size_t point1Index = twoViewGeometry.inlierMatches[i].point2DIndex1;
		const size_t point2Index = twoViewGeometry.inlierMatches[i].point2DIndex2;
		const CKeypoint& keypoint1 = keypoints1[point1Index];
		const CKeypoint& keypoint2 = keypoints2[point2Index];
		const Eigen::Vector3f p1(keypoint1.pt.x, keypoint1.pt.y, 1.0f);
		if (twoViewGeometry.type == CTwoViewGeometryType::CCalibrated || twoViewGeometry.type == CTwoViewGeometryType::CUncalibrated)
		{
			const Eigen::Vector3f p2(keypoint2.pt.x, keypoint2.pt.y, 1.0f);
			const Eigen::Vector3f Fx1 = F * p1;
			const Eigen::Vector3f Ftx2 = F.transpose() * p2;
			const float x2tFx1 = p2.transpose() * Fx1;
			if (x2tFx1 * x2tFx1 / (Fx1(0) * Fx1(0) + Fx1(1) * Fx1(1) + Ftx2(0) * Ftx2(0) + Ftx2(1) * Ftx2(1)) <= max_residual)
			{
				newInlierMatches.push_back(twoViewGeometry.inlierMatches[i]);
			}
		}
		else
		{
			const Eigen::Vector2f p2(keypoint2.pt.x, keypoint2.pt.y);
			if (((H * p1).hnormalized() - p2).squaredNorm() <= max_residual)
			{
				newInlierMatches.push_back(twoViewGeometry.inlierMatches[i]);
			}
		}
	}
	newInlierMatches.shrink_to_fit();
	twoViewGeometry.inlierMatches = newInlierMatches;


	// 原本的策略. 实际测试发现, 使用了guided match之后, 误匹配甚至有可能更多.
	/*
	twoViewGeometry.inlierMatches.clear();
	FlannIndexType flannIndex1 = BuildFlannIndex(descriptors1);
	FlannIndexType flannIndex2 = BuildFlannIndex(descriptors2);
	float max_residual = options.maxError * options.maxError;

	const Eigen::Matrix3f F = twoViewGeometry.F.cast<float>();
	const Eigen::Matrix3f H = twoViewGeometry.H.cast<float>();

	function<bool(float, float, float, float)> guidedFilter;
	if (twoViewGeometry.type == CTwoViewGeometryType::CCalibrated || twoViewGeometry.type == CTwoViewGeometryType::CUncalibrated)
	{
		guidedFilter = [&](const float x1, const float y1, const float x2, const float y2) 
			{
				const Eigen::Vector3f p1(x1, y1, 1.0f);
				const Eigen::Vector3f p2(x2, y2, 1.0f);
				const Eigen::Vector3f Fx1 = F * p1;
				const Eigen::Vector3f Ftx2 = F.transpose() * p2;
				const float x2tFx1 = p2.transpose() * Fx1;
				return x2tFx1 * x2tFx1 / (Fx1(0) * Fx1(0) + Fx1(1) * Fx1(1) + Ftx2(0) * Ftx2(0) + Ftx2(1) * Ftx2(1)) > max_residual;
			};
	}
	else
	{
		guidedFilter = [&](const float x1, const float y1, const float x2, const float y2)
			{
				const Eigen::Vector3f p1(x1, y1, 1.0f);
				const Eigen::Vector2f p2(x2, y2);
				return ((H * p1).hnormalized() - p2).squaredNorm() > max_residual;
			};
	}
	CHECK(guidedFilter != nullptr);
	const Eigen::MatrixXi distances = ComputeSIFTDistanceMatrix(keypoints1, keypoints2, descriptors1, descriptors2, guidedFilter);
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>indices_1to2(distances.rows(), distances.cols());
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>indices_2to1(distances.cols(), distances.rows());
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>distances_1to2 = distances;
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>distances_2to1 = distances.transpose();
	for (int i = 0; i < indices_1to2.rows(); i++)
	{
		indices_1to2.row(i) = Eigen::VectorXi::LinSpaced(indices_1to2.cols(), 0, indices_1to2.cols() - 1);
	}
	for (int i = 0; i < indices_2to1.rows(); i++)
	{
		indices_2to1.row(i) = Eigen::VectorXi::LinSpaced(indices_2to1.cols(), 0, indices_2to1.cols() - 1);
	}
	twoViewGeometry.inlierMatches = FindBestMatchesFlann(indices_1to2, distances_1to2, indices_2to1, distances_2to1, options.maxRatio, options.maxDistance, options.doCrossCheck);
	*/

}
CSIFTCPUMatcher::FlannIndexType CSIFTCPUMatcher::BuildFlannIndex(const CSIFTDescriptors& descriptors)
{
	const flann::Matrix<uint8_t> descriptors_matrix(const_cast<uint8_t*>(descriptors.data()), descriptors.rows(), 128);
	constexpr size_t kNumTreesInForest = 4;

	FlannIndexType index(descriptors_matrix, flann::KDTreeIndexParams(kNumTreesInForest));
	index.buildIndex();
	return index;
}

size_t CSIFTGPUMatcher::lastDescriptors1Index = numeric_limits<size_t>::max();
size_t CSIFTGPUMatcher::lastDescriptors2Index = numeric_limits<size_t>::max();
bool CSIFTGPUMatcher::isUploadDescriptors1 = true;
bool CSIFTGPUMatcher::isUploadDescriptors2 = true;
CSIFTGPUMatcher::CSIFTGPUMatcher(const CSIFTMatchingOptions& options)
{
	options.Check();
	this->options = options;

	SiftGPU siftGPU;
	siftGPU.SetVerbose(0);
	siftMatchGPU = SiftMatchGPU(options.maxNumMatches);

	int gpuIndex = SetBestCudaDevice();
	siftMatchGPU.SetLanguage(SiftMatchGPU::SIFTMATCH_CUDA_DEVICE0 + gpuIndex);
	CHECK(siftMatchGPU.VerifyContextGL() != 0);

	CHECK(siftMatchGPU.Allocate(options.maxNumMatches, options.doCrossCheck), "Allocate ERROR: Not enough GPU memory!");
	siftMatchGPU.gpu_index = gpuIndex;

	lastDescriptors1Index = numeric_limits<size_t>::max();
	lastDescriptors2Index = numeric_limits<size_t>::max();
	isUploadDescriptors1 = true;
	isUploadDescriptors2 = true;
}
CSIFTMatches CSIFTGPUMatcher::Match(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2)
{
	CHECK(descriptors1.cols() == 128 && descriptors2.cols() == 128);
	CHECK(descriptors1.rows() != 0 && descriptors2.rows() != 0);

	struct Uint32_Match
	{
		uint32_t point1Index;
		uint32_t point2Index;
	};
	vector<Uint32_Match> matchesInternal(options.maxNumMatches);
	if (isUploadDescriptors1)
	{
		siftMatchGPU.SetDescriptors(0, descriptors1.rows(), descriptors1.data());
	}
	if (isUploadDescriptors2)
	{
		siftMatchGPU.SetDescriptors(1, descriptors2.rows(), descriptors2.data());
	}
	isUploadDescriptors1 = isUploadDescriptors2 = true; // 每Match一次之后都要求外部程序重新判断是否要上传新的descriptors 0.7f
	const int numMatches = siftMatchGPU.GetSiftMatch(options.maxNumMatches, reinterpret_cast<uint32_t(*)[2]>(matchesInternal.data()), 0.2f, options.maxRatio, options.doCrossCheck);
	if (numMatches < 0)
	{
		cerr << "ERROR: Feature matching failed!" << endl;
		OutputDebugStringA("ERROR: Feature matching failed!\n");
		return CSIFTMatches(0);
	}
	CHECK(numMatches <= matchesInternal.size());
	matchesInternal.resize(numMatches);

	CSIFTMatches matches(numMatches);
	for (size_t i = 0; i < numMatches; i++)
	{
		matches[i].point2DIndex1 = matchesInternal[i].point1Index;
		matches[i].point2DIndex2 = matchesInternal[i].point2Index;
	}
	return matches;
}
void CSIFTGPUMatcher::MatchGuided(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, CTwoViewGeometry& twoViewGeometry)
{
	CHECK(descriptors1.cols() == 128 && descriptors2.cols() == 128);
	CHECK(descriptors1.rows() != 0 && descriptors2.rows() != 0);
	CHECK(keypoints1.size() == descriptors1.rows() && keypoints2.size() == descriptors2.rows());
	if (twoViewGeometry.type < 2 || twoViewGeometry.type > 6)
	{
		return;
	}
	struct Float_Keypoint
	{
		float x;
		float y;

		float a11;
		float a12;
		float a21;
		float a22;
	};
	vector<Float_Keypoint> points1Internal(keypoints1.size()), points2Internal(keypoints2.size());
	for (size_t i = 0; i < keypoints1.size(); i++)
	{
		points1Internal[i].x = keypoints1[i].pt.x;
		points1Internal[i].y = keypoints1[i].pt.y;
		points1Internal[i].a11 = keypoints1[i].a11;
		points1Internal[i].a12 = keypoints1[i].a12;
		points1Internal[i].a21 = keypoints1[i].a21;
		points1Internal[i].a22 = keypoints1[i].a22;
	}
	for (size_t i = 0; i < keypoints2.size(); i++)
	{
		points2Internal[i].x = keypoints2[i].pt.x;
		points2Internal[i].y = keypoints2[i].pt.y;
		points2Internal[i].a11 = keypoints2[i].a11;
		points2Internal[i].a12 = keypoints2[i].a12;
		points2Internal[i].a21 = keypoints2[i].a21;
		points2Internal[i].a22 = keypoints2[i].a22;
	}
	static_assert(offsetof(Float_Keypoint, x) == 0 * sizeof(float), "Invalid keypoint format");
	static_assert(offsetof(Float_Keypoint, y) == 1 * sizeof(float), "Invalid keypoint format");
	static_assert(sizeof(Float_Keypoint) == 6 * sizeof(float), "Invalid keypoint format");

	if (isUploadDescriptors1)
	{
		siftMatchGPU.SetDescriptors(0, descriptors1.rows(), descriptors1.data());
		siftMatchGPU.SetFeautreLocation(0, reinterpret_cast<const float*>(points1Internal.data()), 4);
	}
	if (isUploadDescriptors2)
	{
		siftMatchGPU.SetDescriptors(1, descriptors2.rows(), descriptors2.data());
		siftMatchGPU.SetFeautreLocation(1, reinterpret_cast<const float*>(points2Internal.data()), 4);
	}
	Eigen::Matrix<float, 3, 3, Eigen::RowMajor> F;
	Eigen::Matrix<float, 3, 3, Eigen::RowMajor> H;
	float* F_ptr = nullptr;
	float* H_ptr = nullptr;
	if (twoViewGeometry.type == CTwoViewGeometryType::CCalibrated || twoViewGeometry.type == CTwoViewGeometryType::CUncalibrated)
	{
		F = twoViewGeometry.F.cast<float>();
		F_ptr = F.data();
	}
	else
	{
		H = twoViewGeometry.H.cast<float>();
		H_ptr = H.data();
	}
	CHECK(F_ptr != nullptr || H_ptr != nullptr);

	struct Uint32_Match
	{
		uint32_t point1Index;
		uint32_t point2Index;
	};
	vector<Uint32_Match> matchesInternal(options.maxNumMatches);
	float maxError = options.maxError * options.maxError;
	const int numMatches = siftMatchGPU.GetGuidedSiftMatch(options.maxNumMatches, reinterpret_cast<uint32_t(*)[2]>(matchesInternal.data()), H_ptr, F_ptr, 0.7f, options.maxRatio, maxError, maxError, options.doCrossCheck);
	if (numMatches < 0)
	{
		cerr << "ERROR: Guided feature matching failed!" << endl;
		OutputDebugStringA("ERROR: Guided feature matching failed!\n");
		return;
	}
	CHECK(numMatches <= options.maxNumMatches);
	matchesInternal.resize(numMatches);
	twoViewGeometry.inlierMatches.resize(numMatches);
	for (size_t i = 0; i < numMatches; i++)
	{
		twoViewGeometry.inlierMatches[i].point2DIndex1 = matchesInternal[i].point1Index;
		twoViewGeometry.inlierMatches[i].point2DIndex2 = matchesInternal[i].point2Index;
	}
}


