#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <flann/flann.hpp>
#include "../../src/Scene/Point2D.h"
#include "../../src/Base/Math.h"
#include "../../src/Base/Base.h"

// ��CKeypointsת��Ϊvector<Eigen::Vector2d>
std::vector<Eigen::Vector2d> KeypointsToVector2d(const CKeypoints& keypoints);

// �������ӽ���L2��һ��, ÿһ�б�ʾһ������
void L2Normalize(CSIFTDescriptors_Float& descriptors);

// �������ӽ���L1-Root��һ��, ÿһ�б�ʾһ������. �ο�����: "Three things everyone should know to improve object retrieval", Relja Arandjelovic and Andrew Zisserman, CVPR 2012.
void L1RootNormalize(CSIFTDescriptors_Float& descriptors);

// ��CSIFTDescriptors_Float���͵�������ת����CSIFTDescriptors����, ��[0, 0.5]��Χ���Ա任��[0, 255]. 
// �ضϵ����ֵ0.5��Ϊ�˱��⾫����ʧ, ����ѭ�˱�ʾSIFT�����ĳ�������
CSIFTDescriptors SIFTDescriptorsFloatToUnsignedChar(const CSIFTDescriptors_Float& descriptors);

// ��ȡ�����߶��������Ӧ��������
void ExtractTopScaleFeatures(CKeypoints& keypoints, CSIFTDescriptors& descriptors, size_t numFeatures);

// VLFeatʹ�ò�ͬ��Լ�����洢��������, ���������VLFeat��ʽת��ΪԭʼSIFT��ʽ
CSIFTDescriptors VLFeatToOriginFormat(const CSIFTDescriptors& vlfeatDescriptors);

// ��ȡӰ������, �����ͨ��(�Ҷ�)����vector<uchar>
std::vector<uchar> ReadImage(const std::string& imagePath, size_t maxSize, size_t& originWidth, size_t& originHeight, size_t& currentWidth, size_t& currentHeight);

// ��������ĳ߶ȷֱ���X��Y�����ϷŴ�ԭ����scaleX��scaleY��
void ScaleKeypoints(CKeypoints& keypoints, double scaleX, double scaleY);
void ScaleKeypoints(CKeypoints& keypoints, size_t currentWidth, size_t currentHeight, size_t targetWidth, size_t targetHeight);

// ��������Ӱ�������SIFT������֮��ľ������, ����descriptors1��m��, descriptors2��n��, ��ô���صľ�����m��n��
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2);
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, const std::function<bool(float, float, float, float)>& guidedFilter);

// ���ڸ����ľ�������ҵ���ѵ�����ƥ��
size_t FindBestMatchesOneWayBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, std::vector<int>& matches);

// ���н�����֤��Ѱ���������ƥ��
void FindBestMatchesBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, bool isCrossCheck, CSIFTMatches& matches);

// ʹ��Flann������K��������
void FindNearestNeighborsFlann(const CSIFTDescriptors& query, const CSIFTDescriptors& index, const flann::Index<flann::L2<uint8_t>>& flannIndex, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances);
size_t FindBestMatchesOneWayFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances, float maxRatio, float maxDistance, std::vector<int>& matches);
CSIFTMatches FindBestMatchesFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_2to1, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_2to1, float maxRatio, float maxDistance, bool doCrossCheck);

void ExportMatches(const std::string& leftImagePath, const std::string& rightImagePath, const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTMatches& matches, const std::string& exportImagePath);


