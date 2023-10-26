#include "Common.h"
#include <array>

using namespace std;
vector<Eigen::Vector2d> KeypointsToVector2d(const CKeypoints& keypoints)
{
	vector<Eigen::Vector2d> results(keypoints.size());
	for (size_t i = 0; i < keypoints.size(); i++)
	{
		results[i] = Eigen::Vector2d(keypoints[i].pt.x, keypoints[i].pt.y);
	}
	return results;
}
void L2Normalize(CSIFTDescriptors_Float& descriptors)
{
	descriptors.rowwise().normalize();
}
void L1RootNormalize(CSIFTDescriptors_Float& descriptors)
{
	for (Eigen::MatrixXf::Index r = 0; r < descriptors.rows(); r++)
	{
		descriptors.row(r) *= 1 / descriptors.row(r).lpNorm<1>();
		descriptors.row(r) = descriptors.row(r).array().sqrt();
	}
}
CSIFTDescriptors SIFTDescriptorsFloatToUnsignedChar(const CSIFTDescriptors_Float& descriptors)
{
	CSIFTDescriptors results(descriptors.rows(), descriptors.cols());
	for (Eigen::MatrixXf::Index r = 0; r < descriptors.rows(); r++)
	{
		for (Eigen::MatrixXf::Index c = 0; c < descriptors.cols(); c++)
		{
			float originValue = descriptors(r, c);
			Check(originValue >= 0 && originValue <= 0.5);
			const float scaledValue = round(512.0f * originValue);
			results(r, c) = TruncateCast<float, uint8_t>(scaledValue);
		}
	}
	return results;
}
void ExtractTopScaleFeatures(CKeypoints& keypoints, CSIFTDescriptors& descriptors, size_t numFeatures)
{
	Check(keypoints.size() == descriptors.rows());
	Check(numFeatures > 0);

	if (descriptors.rows() <= numFeatures)
	{
		return;
	}
	vector<pair<size_t, float>> scales;
	scales.reserve(keypoints.size());
	for (size_t i = 0; i < keypoints.size(); i++)
	{
		scales.emplace_back(i, keypoints[i].GetScale());
	}
	partial_sort(scales.begin(), scales.begin() + numFeatures, scales.end(),
		[](const pair<size_t, float>& scale1, const pair<size_t, float>& scale2)
		{
			return scale1.second > scale2.second;
		});
	CKeypoints topScaledKeypoints(numFeatures);
	CSIFTDescriptors topScaleDescriptors(numFeatures, descriptors.cols());
	for (size_t i = 0; i < numFeatures; i++) 
	{
		topScaledKeypoints[i] = keypoints[scales[i].first];
		topScaleDescriptors.row(i) = descriptors.row(scales[i].first);
	}

	keypoints = move(topScaledKeypoints);
	descriptors = move(topScaleDescriptors);
}
CSIFTDescriptors VLFeatToOriginFormat(const CSIFTDescriptors& vlfeatDescriptors)
{
	CSIFTDescriptors originFormat(vlfeatDescriptors.rows(), vlfeatDescriptors.cols());
	const array<int, 8> q{ {0, 7, 6, 5, 4, 3, 2, 1} };
	for (CSIFTDescriptors::Index n = 0; n < vlfeatDescriptors.rows(); n++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 8; k++)
				{
					originFormat(n, 8 * (j + 4 * i) + q[k]) = vlfeatDescriptors(n, 8 * (j + 4 * i) + k);
				}
			}
		}
	}
	return originFormat;
}
vector<uchar> ReadImage(const string& imagePath, size_t maxSize, size_t& originWidth, size_t& originHeight, size_t& currentWidth, size_t& currentHeight)
{
	Check(IsFileExists(imagePath));
	cv::Mat image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
	Check(!image.empty());

	originWidth = image.cols;
	originHeight = image.rows;
	if (image.rows > maxSize || image.cols > maxSize)
	{
		double scale = min(maxSize * 1.0/ image.rows, maxSize * 1.0 / image.cols);
		size_t newRows = image.rows * scale;
		size_t newCols = image.cols * scale;
		cv::resize(image, image, cv::Size(newCols, newRows));
	}

	vector<uchar> imageData;
	imageData.reserve(image.rows * image.cols);
	if (image.isContinuous())
	{
		imageData.assign(image.data, image.data + image.total());
	}
	else
	{
		for (size_t i = 0; i < image.rows; i++)
		{
			imageData.insert(imageData.end(), image.ptr<uchar>(i), image.ptr<uchar>(i) + image.cols);
		}
	}
	currentWidth = image.cols;
	currentHeight = image.rows;
	return imageData;
}
void ScaleKeypoints(CKeypoints& keypoints, double scaleX, double scaleY)
{
	for (size_t i = 0; i < keypoints.size(); i++)
	{
		keypoints[i].Rescale(scaleX, scaleY);
	}
}
void ScaleKeypoints(CKeypoints& keypoints, size_t currentWidth, size_t currentHeight, size_t targetWidth, size_t targetHeight)
{
	double scaleX = targetWidth * 1.0 / currentWidth;
	double scaleY = targetHeight * 1.0 / currentHeight;
	ScaleKeypoints(keypoints, scaleX, scaleY);
}
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2)
{
	Check(descriptors1.cols() == 128 && descriptors2.cols() == 128);
	Check(descriptors1.rows() != 0 && descriptors2.rows() != 0);
	int m = descriptors1.rows(), n = descriptors2.rows();

	Eigen::MatrixXi distances(m, n);
#pragma omp parallel for
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			distances(i, j) = (descriptors1.row(i).cast<float>() - descriptors2.row(j).cast<float>()).norm();
		}
	}
	return distances;
}
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, const function<bool(float, float, float, float)>& guidedFilter)
{
	Check(guidedFilter != nullptr);
	Check(keypoints1.size() == descriptors1.rows() && keypoints2.size() == descriptors2.rows());
	int m = descriptors1.rows(), n = descriptors2.rows();
	Eigen::MatrixXi distances(m, n);
#pragma omp parallel for
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (guidedFilter(keypoints1[i].pt.x, keypoints1[i].pt.y, keypoints2[j].pt.x, keypoints2[j].pt.y))
			{
				distances(i, j) = numeric_limits<int>::max();
			}
			else
			{
				distances(i, j) = (descriptors1.row(i).cast<float>() - descriptors2.row(j).cast<float>()).norm();
			}
		}
	}
	return distances;
}
size_t FindBestMatchesOneWayBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, vector<int>& matches)
{
	size_t numMatches = 0;
	matches.resize(distances.rows(), -1);

	for (size_t index1 = 0; index1 < distances.rows(); index1++)
	{
		int bestIndex2 = -1;
		int minDistance = numeric_limits<int>::max();
		int secondMinDistance = numeric_limits<int>::max();

		for (size_t index2 = 0; index2 < distances.cols(); index2++)
		{
			const int dist = distances(index1, index2);
			if (dist < minDistance)
			{
				bestIndex2 = index2;
				secondMinDistance = minDistance;
				minDistance = dist;
			}
			else if (dist < secondMinDistance) 
			{
				secondMinDistance = dist;
			}
		}
		if (bestIndex2 == -1) 
		{
			continue;
		}
		if (minDistance > maxDistance)
		{
			continue;
		}
		if (minDistance >= maxRatio * secondMinDistance)
		{
			continue;
		}
		numMatches++;
		matches[index1] = bestIndex2;
	}
	return numMatches;
}
void FindBestMatchesBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, bool isCrossCheck, CSIFTMatches& matches)
{
	matches.clear();
	vector<int> matches12;
	const size_t numMatches12 = FindBestMatchesOneWayBruteForce(distances, maxRatio, maxDistance, matches12);
	if (isCrossCheck)
	{
		vector<int> matches21;
		const size_t numMatches21 = FindBestMatchesOneWayBruteForce(distances.transpose(), maxRatio, maxDistance, matches21);
		for (size_t i1 = 0; i1 < matches12.size(); i1++)
		{
			if (matches12[i1] != -1 && matches21[matches12[i1]] != -1 && matches21[matches12[i1]] == i1)
			{
				matches.push_back(CSIFTMatch(i1, matches12[i1]));
			}
		}
	}
	else
	{
		matches.reserve(numMatches12);
		for (size_t i1 = 0; i1 < matches12.size(); i1++)
		{
			if (matches12[i1] != -1)
			{
				matches.push_back(CSIFTMatch(i1, matches12[i1]));
			}
		}
	}
}
void FindNearestNeighborsFlann(const CSIFTDescriptors& query, const CSIFTDescriptors& index, const flann::Index<flann::L2<uint8_t>>& flannIndex, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances)
{
	if (query.rows() == 0 || index.rows() == 0)
	{
		return;
	}

	size_t numNearestNeighbors = 2;           // 要查找的最近邻的数量
	constexpr size_t numLeafsToVisit = 128;   // FLANN搜索时要访问的叶节点的最大数量
	numNearestNeighbors = min(numNearestNeighbors, static_cast<size_t>(index.rows()));

	indices.resize(query.rows(), numNearestNeighbors);
	distances.resize(query.rows(), numNearestNeighbors);
	const flann::Matrix<uint8_t> queryMatrix(const_cast<uint8_t*>(query.data()), query.rows(), 128);

	flann::Matrix<int> indicesMatrix(indices.data(), query.rows(), numNearestNeighbors);
	vector<float> distancesVector(query.rows() * numNearestNeighbors);
	flann::Matrix<float> distancesMatrix(distancesVector.data(), query.rows(), numNearestNeighbors);
	flannIndex.knnSearch(queryMatrix, indicesMatrix, distancesMatrix, numNearestNeighbors, flann::SearchParams(numLeafsToVisit));

	for (Eigen::Index queryIndex = 0; queryIndex < indices.rows(); queryIndex++)
	{
		for (Eigen::Index k = 0; k < indices.cols(); k++)
		{
			const Eigen::Index indexIndex = indices.coeff(queryIndex, k);
			distances(queryIndex, k) = (query.row(queryIndex).cast<float>() - index.row(indexIndex).cast<float>()).norm();
		}
	}
}
size_t FindBestMatchesOneWayFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances, float maxRatio, float maxDistance, std::vector<int>& matches)
{
	size_t numMatches = 0;
	matches.resize(indices.rows(), -1);
	for (int index1 = 0; index1 < indices.rows(); index1++)
	{
		int bestIndex2 = -1;
		int minDistance = numeric_limits<int>::max();
		int secondMinDistance = numeric_limits<int>::max();
		for (int nIndex = 0; nIndex < indices.cols(); nIndex++)
		{
			const int index2 = indices(index1, nIndex);
			const int dist = distances(index1, nIndex);
			if (dist < minDistance)
			{
				bestIndex2 = index2;
				secondMinDistance = minDistance;
				minDistance = dist;
			}
			else if (dist < secondMinDistance)
			{
				secondMinDistance = dist;
			}
		}
		if (bestIndex2 == -1)
		{
			continue;
		}
		if (minDistance > maxDistance)
		{
			continue;
		}
		if (minDistance >= maxRatio * secondMinDistance)
		{
			continue;
		}
		numMatches++;
		matches[index1] = bestIndex2;
	}
	return numMatches;
}
CSIFTMatches FindBestMatchesFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_2to1, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_2to1, float maxRatio, float maxDistance, bool doCrossCheck)
{
	CSIFTMatches matches;
	vector<int> matches12;
	const size_t numMatches12 = FindBestMatchesOneWayFlann(indices_1to2, distances_1to2, maxRatio, maxDistance, matches12);
	if (doCrossCheck && indices_2to1.rows())
	{
		vector<int> matches21;
		const size_t numMatches21 = FindBestMatchesOneWayFlann(indices_2to1, distances_2to1, maxRatio, maxDistance, matches21);
		matches.reserve(min(numMatches12, numMatches21));
		for (size_t i1 = 0; i1 < matches12.size(); i1++)
		{
			if (matches12[i1] != -1 && matches21[matches12[i1]] != -1 && matches21[matches12[i1]] == static_cast<int>(i1))
			{
				matches.push_back(CSIFTMatch(i1, matches12[i1]));
			}
		}
	}
	else
	{
		matches.reserve(numMatches12);
		for (size_t i1 = 0; i1 < matches12.size(); i1++)
		{
			if (matches12[i1] != -1)
			{
				matches.push_back(CSIFTMatch(i1, matches12[i1]));
			}
		}
	}
	return matches;
}

void ExportMatches(const string& leftImagePath, const string& rightImagePath, const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTMatches& matches, const string& exportImagePath)
{
	const cv::Mat leftImage = cv::imread(leftImagePath, cv::IMREAD_COLOR);
	const cv::Mat rightImage = cv::imread(rightImagePath, cv::IMREAD_COLOR);

	vector<cv::KeyPoint> keypointsCV1(keypoints1.size());
	for (size_t i = 0; i < keypoints1.size(); i++)
	{
		keypointsCV1[i] = keypoints1[i];
	}
	vector<cv::KeyPoint> keypointsCV2(keypoints2.size());
	for (size_t i = 0; i < keypoints2.size(); i++)
	{
		keypointsCV2[i] = keypoints2[i];
	}

	cv::Mat outputImage;
	const cv::Size size(leftImage.cols + rightImage.cols, std::max(leftImage.rows, rightImage.rows));
	outputImage.create(size, CV_MAKETYPE(leftImage.depth(), 3));
	outputImage = cv::Scalar::all(0);

	leftImage.copyTo(outputImage(cv::Rect(0, 0, leftImage.cols, leftImage.rows)));
	rightImage.copyTo(outputImage(cv::Rect(leftImage.cols, 0, rightImage.cols, rightImage.rows)));

	for (const auto& match : matches)
	{
		const cv::Point2f pt1 = keypoints1[match.point2DIndex1].pt;
		cv::Point2f pt2 = keypoints2[match.point2DIndex2].pt;
		pt2.x += leftImage.cols;
		cv::line(outputImage, pt1, pt2, cv::Scalar(0, 255, 0), 4);
	}

	cv::imwrite(exportImagePath, outputImage);
}







