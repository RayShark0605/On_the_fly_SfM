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

	// ׼���ؽ���ģ��
	void BeginReconstruction(const COptions& options);

	// �����ؽ���ģ�Ͳ���ָʾģ���Ƿ�Ӧ�ñ�����
	void EndReconstruction(bool isDiscard);

	// Ѱ�ҳ�ʼӰ���
	bool FindInitialImagePair(const COptions& options, size_t& imageID1, size_t& imageID2);

	// ע���ʼӰ���
	bool RegisterInitialImagePair(const COptions& options, size_t imageID1, size_t imageID2);

	// Ѱ����ѵ���һ��Ӱ��ȥע��
	std::vector<size_t> FindNextImages(const COptions& options);

	// ע����һ��Ӱ��
	bool RegisterNextImage(const COptions& options, size_t imageID);

	// ��Ӱ��Ĺ۲����ǰ������
	size_t TriangulateImage(const COptions& options, size_t imageID);

	// �������ǲ�����Щ������Ӧ��ϵӦ���й�ͬ�۲�㵫����Ư�Ƶ�ԭ���û�е�Ӱ���.
	// Ϊ�˴���Ư��, ��ʹ�õ���ͶӰ�����ֵӦ����Խϴ�. ��ֵ���󽫻ᵼ�·�³����ƽ��ʧЧ, ��ֵ��С���޷���Ч������Ư��.
	size_t Retriangulate(const COptions& options);

	// ���ݵظ�����Ӧ��ϵ����ɹ켣. ��ƽ��֮��, ��һ���ǳ�����, ��Ϊ��������3D��λ���Ѿ����Ż�, ��ɹ켣�ܸ��õ�֧�ֺ�����Ӱ���ע��
	size_t CompleteTracks(const COptions& options);

	// ͨ����Ӧ��ϵ���ϲ��켣. ��CompleteTracks����, ����ƽ��֮��ǳ���Ч, �ܸ��ƺ���ƽ���Ż��е������
	size_t MergeTracks(const COptions& options);

	// �Ż���imageID�ֲ����ӵ�Ӱ��͵�. ����, �Ż������3D��. ֻ����imageID���ӵ�Ӱ��Żᱻ�Ż�.
	// ��������3D����imageIDû�оֲ�����, ��ô�۲���Щ���Ӱ����ƽ������лᱻ����Ϊ����
	CLocalBundleAdjustmentReport LocalBundleAdjust(const COptions& options, size_t imageID, const std::unordered_set<size_t>& point3DIDs);

	// ȫ��ƽ��
	bool GlobalBundleAdjust(const COptions& options);

	// ����Ӱ��͵�۲�ֵ
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

	// ���һ�ε���FindFirstInitialImageʱ���Ƶ�˫�Ӽ���, ����֮���RegisterInitialImagePairʱ�Ļ���
	std::pair<size_t, size_t> preInitImagePair;
	CTwoViewGeometry preInitTwoViewGeometry;

	// �ѳ��Գ�ʼ����Ӱ���, ÿ��Ӱ���ֻ����һ�ν��г�ʼ��
	std::unordered_set<std::pair<size_t, size_t>, MatchPairHash, MatchPairEqual> initImagePairs;

	// ÿ���������ע��ͼ������. ���ڱ����ڱ��ؾֲ�ƽ���ж��Ӱ�����ڲ�ʱ, ����ڲε��ظ�ϸ������ϸ������������˻�
	std::unordered_map<size_t, size_t> numRegImagesPerCamera;

	// ��ǰģ���б����˵�Ӱ��
	std::unordered_set<size_t> filteredImages;

	// �ڿ�ʼ�ؽ�֮ǰ���Ѿ�ע���Ӱ��, һ���Ǽ���������ģ��֮��ŷǿ�
	std::unordered_set<size_t> existingImages;

	// Ѱ�ҳ�ʼӰ��Եĵ�һ��Ӱ��. ����ֵ�ǰ��պ��ʶ������, ����ʵ�Ӱ��λ����ǰ��
	std::vector<size_t> FindFirstInitialImage() const;

	// ���ڸ����ĵ�һ��Ӱ��, Ѱ�ҳ�ʼӰ��Եĵڶ���Ӱ��. ���ʵĵڶ���Ӱ��Ӧ�����һ��Ӱ���д�����ƥ���ϵ, ����������������Ϣ.
	// ����ֵ�ǰ��պ��ʶ������, ����ʵ�Ӱ��λ����ǰ��
	std::vector<size_t> FindSecondInitialImage(const COptions& options, size_t imageID1) const;

	// ���Ҳ���ֲ�ƽ���Ӱ��
	std::vector<size_t> FindLocalBundle(const COptions& options, size_t imageID) const;

	// �ڵ�ǰģ����ע���ȡ��ע��Ӱ��
	void RegisterImageEvent(size_t imageID);
	void DeRegisterImageEvent(size_t imageID);

	// ���Ƴ�ʼӰ��Ե�˫�Ӽ��ι�ϵ
	bool EstimateInitialTwoViewGeometry(const COptions& options, size_t imageID1, size_t imageID2);
};


























