#pragma once
#include "../Scene/Model.h"
#include "../Geometry/Triangulation.h"

struct CLocalBundleAdjustmentReport
{
	size_t numMergedObservations = 0;
	size_t numCompletedObservations = 0;
	size_t numFilteredObservations = 0;
	size_t numAdjustedObservations = 0;
};
class CMapper final
{
public:
	CMapper(CModel& model, CDatabase* const database);
	~CMapper();

	// 准备重建该模型
	void BeginReconstruction(const COptions& options);

	// 结束重建该模型并且指示模型是否应该被丢弃
	void EndReconstruction(bool isDiscard);

	// 寻找初始影像对
	bool FindInitialImagePair(const COptions& options, size_t& imageID1, size_t& imageID2);

	// 注册初始影像对
	bool RegisterInitialImagePair(const COptions& options, size_t imageID1, size_t imageID2);

	// 寻找最佳的下一张影像去注册
	std::vector<size_t> FindNextImages(const COptions& options);

	// 注册下一张影像
	bool RegisterNextImage(const COptions& options, size_t imageID);

	// 对影像的观测点做前方交会
	size_t TriangulateImage(const COptions& options, size_t imageID);

	// 重新三角测量那些根据响应关系应该有共同观测点但由于漂移等原因而没有的影像对.
	// 为了处理漂移, 所使用的重投影误差阈值应该相对较大. 阈值过大将会导致非鲁棒的平差失效, 阈值过小则无法有效地修正漂移.
	size_t Retriangulate(const COptions& options);

	// 传递地根据响应关系来完成轨迹. 在平差之后, 这一步非常有用, 因为许多相机和3D点位置已经被优化, 完成轨迹能更好地支持后续新影像的注册
	size_t CompleteTracks(const COptions& options);

	// 通过响应关系来合并轨迹. 与CompleteTracks类似, 这在平差之后非常有效, 能改善后续平差优化中的冗余度
	size_t MergeTracks(const COptions& options);

	// 优化与imageID局部连接的影像和点. 此外, 优化传入的3D点. 只有与imageID连接的影像才会被优化.
	// 如果传入的3D点与imageID没有局部连接, 那么观测这些点的影像在平差过程中会被设置为常量
	CLocalBundleAdjustmentReport LocalBundleAdjust(const COptions& options, size_t imageID, const std::unordered_set<size_t>& point3DIDs);

	// 全局平差
	bool GlobalBundleAdjust(const COptions& options);

	// 过滤影像和点观测值
	size_t FilterImages(const COptions& options);
	size_t FilterPoints(const COptions& options);

	const CModel& GetModel() const;

	const std::unordered_set<size_t>& GetModifiedPoints3D();

	void ClearModifiedPoints3D();

private:
	CModel& model;
	CTriangulator* const triangulator;
	CDatabase* const database;
	const size_t modelID;

	// 最后一次调用FindFirstInitialImage时估计的双视几何, 用作之后的RegisterInitialImagePair时的缓存
	std::pair<size_t, size_t> preInitImagePair;
	CTwoViewGeometry preInitTwoViewGeometry;

	// 已尝试初始化的影像对, 每个影像对只尝试一次进行初始化
	std::unordered_set<std::pair<size_t, size_t>, MatchPairHash, MatchPairEqual> initImagePairs;

	// 每个相机的已注册图像数量. 用于避免在本地局部平差中多个影像共享内参时, 相机内参的重复细化和已细化相机参数的退化
	std::unordered_map<size_t, size_t> numRegImagesPerCamera;

	// 当前模型中被过滤的影像
	std::unordered_set<size_t> filteredImages;

	// 在开始重建之前就已经注册的影像, 一般是加载了现有模型之后才非空
	std::unordered_set<size_t> existingImages;

	// 寻找初始影像对的第一张影像. 返回值是按照合适度排序的, 最合适的影像位于最前面
	std::vector<size_t> FindFirstInitialImage() const;

	// 对于给定的第一张影像, 寻找初始影像对的第二张影像. 合适的第二张影像应该与第一张影像有大量的匹配关系, 而且有先验的相机信息.
	// 返回值是按照合适度排序的, 最合适的影像位于最前面
	std::vector<size_t> FindSecondInitialImage(const COptions& options, size_t imageID1) const;

	// 查找参与局部平差的影像
	std::vector<size_t> FindLocalBundle(const COptions& options, size_t imageID) const;

	// 在当前模型中注册或取消注册影像
	void RegisterImageEvent(size_t imageID);
	void DeRegisterImageEvent(size_t imageID);

	// 估计初始影像对的双视几何关系
	bool EstimateInitialTwoViewGeometry(const COptions& options, size_t imageID1, size_t imageID2);
};


























