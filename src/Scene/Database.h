#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <tbb/concurrent_unordered_map.h>


#include "../Base/Base.h"
#include "Camera.h"
#include "Image.h"
#include "TwoViewGeometry.h"

class CDatabase
{
public:
	virtual ~CDatabase() {};
	virtual void Clear() = 0;

	virtual size_t GetCamerasNum() const = 0;
	virtual bool IsCameraExists(size_t cameraID) const = 0;
	virtual bool IsCameraExists(const CCamera& camera) const = 0;
	virtual size_t GetCameraID(const CCamera& camera) const = 0;
	virtual const CCamera& GetCamera(size_t cameraID) const = 0;
	virtual CCamera& GetCamera(size_t cameraID) = 0;
	virtual std::vector<CCamera> GetAllCameras() const = 0;
	virtual size_t AddCamera(const CCamera& camera) = 0;

	virtual size_t GetImagesNum() const = 0;
	virtual bool IsImageExists(size_t imageID) const = 0;
	virtual bool IsImageExists(const CImage& image) const = 0;
	virtual bool IsImageExists(const std::string& imageName) const = 0;
	virtual size_t GetImageID(const std::string& imageName) const = 0;
	virtual const CCamera& GetImageCamera(size_t imageID) const = 0;
	virtual const CImage& GetImage(size_t imageID) const = 0;
	virtual CImage& GetImage(size_t imageID) = 0;
	virtual std::vector<CImage> GetAllImages() const = 0;
	virtual std::vector<std::string> GetAllImageNames() const = 0;
	virtual size_t AddImage(const CImage& image) = 0;

	virtual size_t GetImageKeypointsNum(size_t imageID) const = 0;
	virtual const CKeypoints& GetImageKeypoints(size_t imageID) const = 0;
	virtual CKeypoints& GetImageKeypoints(size_t imageID) = 0;
	virtual void AddImageKeypoints(size_t imageID, const CKeypoints& keypoints) = 0;

	virtual size_t GetImageDescriptorsNum(size_t imageID) const = 0;
	virtual const CSIFTDescriptors& GetImageDescriptors(size_t imageID) const = 0;
	virtual CSIFTDescriptors& GetImageDescriptors(size_t imageID) = 0;
	virtual void AddImageDescriptors(size_t imageID, const CSIFTDescriptors& descriptors) = 0;

	virtual size_t GetMatchesNum(size_t imageID1, size_t imageID2) const = 0;
	virtual const CSIFTMatches& GetMatches(size_t imageID1, size_t imageID2) const = 0;
	virtual CSIFTMatches& GetMatches(size_t imageID1, size_t imageID2) = 0;
	virtual std::unordered_map<size_t, std::unordered_map<size_t, CSIFTMatches>> GetAllMatches() const = 0;
	virtual void AddMatches(size_t imageID1, size_t imageID2, const CSIFTMatches& matches) = 0;
	virtual size_t GetMatchedImagesNum(size_t imageID) const = 0;

	virtual size_t GetInlierMatchesNum(size_t imageID1, size_t imageID2) const = 0;
	virtual const CTwoViewGeometry& GetTwoViewGeometry(size_t imageID1, size_t imageID2) const = 0;
	virtual CTwoViewGeometry& GetTwoViewGeometry(size_t imageID1, size_t imageID2) = 0;
	virtual std::unordered_map<size_t, std::unordered_map<size_t, CTwoViewGeometry>> GetAllTwoViewGeometries() const = 0;
	virtual void AddTwoViewGeometry(size_t imageID1, size_t imageID2, const CTwoViewGeometry& twoViewGeometry) = 0;
	virtual size_t GetTwoViewGeometryImagesNum(size_t imageID) const = 0;

	virtual std::unordered_map<std::pair<size_t, size_t>, size_t, MatchPairHash, MatchPairEqual> GetAllCorrespondences() const = 0;

	virtual void SaveAsDir(const std::string& dirPath) const = 0;
};


class CRAMDatabase final :public CDatabase
{
public:
	CRAMDatabase() noexcept;
	~CRAMDatabase();
	void Clear() noexcept override;

	size_t GetCamerasNum() const noexcept override;
	bool IsCameraExists(size_t cameraID) const override;
	bool IsCameraExists(const CCamera& camera) const override;
	size_t GetCameraID(const CCamera& camera) const override;
	const CCamera& GetCamera(size_t cameraID) const override;
	CCamera& GetCamera(size_t cameraID) override;
	std::vector<CCamera> GetAllCameras() const override;
	size_t AddCamera(const CCamera& camera) override;

	size_t GetImagesNum() const noexcept override;
	bool IsImageExists(size_t imageID) const override;
	bool IsImageExists(const CImage& image) const override;
	bool IsImageExists(const std::string& imageName) const override;
	size_t GetImageID(const std::string& imageName) const override;
	const CCamera& GetImageCamera(size_t imageID) const override;
	const CImage& GetImage(size_t imageID) const override;
	CImage& GetImage(size_t imageID) override;
	std::vector<CImage> GetAllImages() const override;
	std::vector<std::string> GetAllImageNames() const override;
	size_t AddImage(const CImage& image) override;

	size_t GetImageKeypointsNum(size_t imageID) const override;
	const CKeypoints& GetImageKeypoints(size_t imageID) const override;
	CKeypoints& GetImageKeypoints(size_t imageID) override;
	void AddImageKeypoints(size_t imageID, const CKeypoints& keypoints) override;

	size_t GetImageDescriptorsNum(size_t imageID) const override;
	const CSIFTDescriptors& GetImageDescriptors(size_t imageID) const override;
	CSIFTDescriptors& GetImageDescriptors(size_t imageID) override;
	void AddImageDescriptors(size_t imageID, const CSIFTDescriptors& descriptors) override;

	size_t GetMatchesNum(size_t imageID1, size_t imageID2) const override;
	const CSIFTMatches& GetMatches(size_t imageID1, size_t imageID2) const override;
	CSIFTMatches& GetMatches(size_t imageID1, size_t imageID2) override;
	std::unordered_map<size_t, std::unordered_map<size_t, CSIFTMatches>> GetAllMatches() const override;
	void AddMatches(size_t imageID1, size_t imageID2, const CSIFTMatches& matches) override;
	size_t GetMatchedImagesNum(size_t imageID) const override;

	size_t GetInlierMatchesNum(size_t imageID1, size_t imageID2) const override;
	const CTwoViewGeometry& GetTwoViewGeometry(size_t imageID1, size_t imageID2) const override;
	CTwoViewGeometry& GetTwoViewGeometry(size_t imageID1, size_t imageID2) override;
	std::unordered_map<size_t, std::unordered_map<size_t, CTwoViewGeometry>> GetAllTwoViewGeometries() const override;
	void AddTwoViewGeometry(size_t imageID1, size_t imageID2, const CTwoViewGeometry& twoViewGeometry) override;
	size_t GetTwoViewGeometryImagesNum(size_t imageID) const override;

	std::unordered_map<std::pair<size_t, size_t>, size_t, MatchPairHash, MatchPairEqual> GetAllCorrespondences() const override;

	void SaveAsDir(const std::string& dirPath) const override;

private:
	tbb::concurrent_unordered_map<size_t, CCamera> cameras;               //key:相机ID, value:相机
	tbb::concurrent_unordered_map<size_t, CImage> images;                 //key:影像ID, value:影像
	//tbb::concurrent_unordered_map<size_t, CKeypoints> keypoints;          //key:影像ID, value:该影像的所有特征点
	//tbb::concurrent_unordered_map<size_t, CSIFTDescriptors> descriptors;  //key:影像ID, value:该影像的所有特征点对应的描述子

	// key:影像ID(较小的那个), value.key:影像ID(较大的那个), value.value:两张影像之间的匹配关系
	tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_map<size_t, CSIFTMatches>> matches;

	// key:影像ID(较小的那个), value.key:影像ID(较大的那个), value.value:两张影像之间的双视几何关系
	tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_map<size_t, CTwoViewGeometry>> twoViewGeometries;

	// 并查集
	class DisjointSet
	{
	public:
		// 查找元素x的根节点, 使用路径压缩的优化机制
		std::pair<size_t, size_t> Find(const std::pair<size_t, size_t>& x)
		{
			if (parent.find(x) == parent.end()) // 如果x不在parent中, 将其作为一个独立集合加入
			{
				parent[x] = x; // x的父节点是其自身
				rank[x] = 0;   // x的秩初始化为0
			}

			// 路径压缩: 如果x不是自己的父节点, 递归地找到其根节点并更新父节点
			if (parent[x] != x)
			{
				parent[x] = Find(parent[x]); // 路径压缩
			}
			return parent[x]; // 返回x的根节点
		}
		// 合并元素x和y所在的集合
		void Unite(std::pair<size_t, size_t> x, std::pair<size_t, size_t> y)
		{
			std::pair<size_t, size_t> xRoot = Find(x); // 找到x的根节点
			std::pair<size_t, size_t> yRoot = Find(y); // 找到y的根节点

			if (xRoot != yRoot) // 如果x和y不在同一个集合中, 则需要合并
			{
				if (rank[xRoot] < rank[yRoot]) // 按秩合并: 将秩较小的根节点的父节点指向秩较大的根节点
				{
					std::swap(xRoot, yRoot);
				}
				parent[yRoot] = xRoot; // y的根节点现在指向x的根节点
				if (rank[xRoot] == rank[yRoot])
				{
					rank[xRoot]++; // 如果秩相等, 则新的根节点的秩加1
				}
			}
		}
		// 获取所有具有相同根节点的元素集合
		std::vector<std::pair<size_t, size_t>> GetConnectedComponents(const std::pair<size_t, size_t>& x)
		{
			std::vector<std::pair<size_t, size_t>> connectedComponents;
			for (const auto& p : parent) // 遍历所有元素, 将与x属于同一集合的元素加入列表
			{
				if (Find(p.first) == Find(x)) // 与x同根的元素
				{
					connectedComponents.push_back(p.first);
				}
			}
			return connectedComponents;
		}

	private:	
		std::unordered_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>, MatchPairHash> parent; // 存储每个元素的父节点
		std::unordered_map<std::pair<size_t, size_t>, int, MatchPairHash> rank;                         // 存储每个根节点的秩
	};
	DisjointSet disjoinSet;
};


















