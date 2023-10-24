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

// 将CKeypoints转换为vector<Eigen::Vector2d>
std::vector<Eigen::Vector2d> KeypointsToVector2d(const CKeypoints& keypoints);

// 对描述子进行L2归一化, 每一行表示一个特征
void L2Normalize(CSIFTDescriptors_Float& descriptors);

// 对描述子进行L1-Root归一化, 每一行表示一个特征. 参考文献: "Three things everyone should know to improve object retrieval", Relja Arandjelovic and Andrew Zisserman, CVPR 2012.
void L1RootNormalize(CSIFTDescriptors_Float& descriptors);

// 将CSIFTDescriptors_Float类型的描述子转换成CSIFTDescriptors类型, 从[0, 0.5]范围线性变换到[0, 255]. 
// 截断到最大值0.5是为了避免精度损失, 并遵循了表示SIFT向量的常见做法
CSIFTDescriptors SIFTDescriptorsFloatToUnsignedChar(const CSIFTDescriptors_Float& descriptors);

// 提取与最大尺度特征相对应的描述子
void ExtractTopScaleFeatures(CKeypoints& keypoints, CSIFTDescriptors& descriptors, size_t numFeatures);

// VLFeat使用不同的约定来存储其描述符, 这个方法将VLFeat格式转换为原始SIFT格式
CSIFTDescriptors VLFeatToOriginFormat(const CSIFTDescriptors& vlfeatDescriptors);

// 读取影像数据, 输出单通道(灰度)数据vector<uchar>
std::vector<uchar> ReadImage(const std::string& imagePath, size_t maxSize, size_t& originWidth, size_t& originHeight, size_t& currentWidth, size_t& currentHeight);

// 将特征点的尺度分别在X和Y方向上放大到原来的scaleX和scaleY倍
void ScaleKeypoints(CKeypoints& keypoints, double scaleX, double scaleY);
void ScaleKeypoints(CKeypoints& keypoints, size_t currentWidth, size_t currentHeight, size_t targetWidth, size_t targetHeight);

// 计算两张影像的所有SIFT特征点之间的距离矩阵, 假设descriptors1有m行, descriptors2有n行, 那么返回的矩阵是m行n列
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2);
Eigen::MatrixXi ComputeSIFTDistanceMatrix(const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTDescriptors& descriptors1, const CSIFTDescriptors& descriptors2, const std::function<bool(float, float, float, float)>& guidedFilter);

// 基于给定的距离矩阵找到最佳的特征匹配
size_t FindBestMatchesOneWayBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, std::vector<int>& matches);

// 带有交叉验证的寻找最佳特征匹配
void FindBestMatchesBruteForce(const Eigen::MatrixXi& distances, float maxRatio, float maxDistance, bool isCrossCheck, CSIFTMatches& matches);

// 使用Flann来加速K近邻搜索
void FindNearestNeighborsFlann(const CSIFTDescriptors& query, const CSIFTDescriptors& index, const flann::Index<flann::L2<uint8_t>>& flannIndex, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances);
size_t FindBestMatchesOneWayFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances, float maxRatio, float maxDistance, std::vector<int>& matches);
CSIFTMatches FindBestMatchesFlann(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_1to2, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices_2to1, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& distances_2to1, float maxRatio, float maxDistance, bool doCrossCheck);

void ExportMatches(const std::string& leftImagePath, const std::string& rightImagePath, const CKeypoints& keypoints1, const CKeypoints& keypoints2, const CSIFTMatches& matches, const std::string& exportImagePath);


