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

	void SaveAsDir(const std::string& dirPath) const override;

private:
	tbb::concurrent_unordered_map<size_t, CCamera> cameras;               //key:���ID, value:���
	tbb::concurrent_unordered_map<size_t, CImage> images;                 //key:Ӱ��ID, value:Ӱ��
	//tbb::concurrent_unordered_map<size_t, CKeypoints> keypoints;          //key:Ӱ��ID, value:��Ӱ�������������
	//tbb::concurrent_unordered_map<size_t, CSIFTDescriptors> descriptors;  //key:Ӱ��ID, value:��Ӱ��������������Ӧ��������

	// key:Ӱ��ID(��С���Ǹ�), value.key:Ӱ��ID(�ϴ���Ǹ�), value.value:����Ӱ��֮���ƥ���ϵ
	tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_map<size_t, CSIFTMatches>> matches;

	// key:Ӱ��ID(��С���Ǹ�), value.key:Ӱ��ID(�ϴ���Ǹ�), value.value:����Ӱ��֮���˫�Ӽ��ι�ϵ
	tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_map<size_t, CTwoViewGeometry>> twoViewGeometries;
};

















