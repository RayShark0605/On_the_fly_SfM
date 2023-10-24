#pragma once
#include <cstddef>
#include <vector>
#include <limits>

#include "../Base/Base.h"
#include "../Base/Math.h"
//#include "Estimator.h"

// ������. ����һ��������
class CSampler
{
public:
	CSampler() = default;
	explicit CSampler(size_t numSamples);
	virtual ~CSampler() = default;

	// ��ʼ��������. �ڵ���Sample֮ǰ��������������.
	virtual void Initialize(size_t numTotalSamples) = 0;

	// ��������Ψһ�������������
	virtual size_t GetMaxNumSamples() const = 0;

	// ���������в���numSamples��Ԫ��, ������������Ԫ�ص������洢��sampledIndex��
	virtual void Sample(std::vector<size_t>& sampledIndex) = 0;

	// ��һϵ��X�в����õ�X_rand
	template <typename XType>
	void GetSampleX(const std::vector<XType>& X, std::vector<XType>& X_rand)
	{
		CHECK(X.size() == numTotalSamples);
		CHECK(numSamples <= numTotalSamples);

		thread_local std::vector<size_t> sampledIndex;
		Sample(sampledIndex);
		X_rand.resize(numSamples);
		for (size_t i = 0; i < numSamples; i++)
		{
			X_rand[i] = X[sampledIndex[i]];
		}
	}

	// ��һϵ��X�в����õ�X_rand, ��һϵ��Y�в����õ�Y_rand
	template <typename XType, typename YType>
	void GetSampleXY(const std::vector<XType>& X, const std::vector<YType>& Y, std::vector<XType>& X_rand, std::vector<YType>& Y_rand)
	{
		CHECK(X.size() == numTotalSamples);
		CHECK(Y.size() == numTotalSamples);

		thread_local std::vector<size_t> sampledIndex;
		Sample(sampledIndex);
		X_rand.resize(numSamples);
		Y_rand.resize(numSamples);
		for (size_t i = 0; i < numSamples; i++)
		{
			X_rand[i] = X[sampledIndex[i]];
			Y_rand[i] = Y[sampledIndex[i]];
		}
	}

protected:
	const size_t numSamples;
	size_t numTotalSamples;
};

// ���������
class CRandomSampler :public CSampler
{
public:
	explicit CRandomSampler(size_t numSamples);
	//explicit CRandomSampler(const CEstimator& estimator);

	void Initialize(size_t numTotalSamples) override;
	size_t GetMaxNumSamples() const override;
	void Sample(std::vector<size_t>& sampledIndex) override;

private:
	std::vector<size_t> allIndex;
};

// PROSAC����(Progressive Sample Consensus)������. �ο�����: "Matching with PROSAC - Progressive Sample Consensus". Ondrej Chum and Matas, CVPR 2005.
// ÿ���̶߳�Ӧ��ʵ����һ�������Ĳ�����, �����������ߵ�����Ӧ�ÿ����б��ǰ��
class CProgressiveSampler : public CSampler
{
public:
	explicit CProgressiveSampler(size_t numSamples);
	//explicit CProgressiveSampler(const CEstimator& estimator);

	void Initialize(size_t numTotalSamples) override;
	size_t GetMaxNumSamples() const override;
	void Sample(std::vector<size_t>& sampledIndex) override;

private:
	size_t t, n;
	double T_n, T_n_p; // ��PROSAC�㷨�ж�����뷽��3��صı���
};

// ���ڻ���RANSAC���������������, �ò���������Ψһ������
// ÿ���߳�Ӧ��ʵ����һ�������Ĳ�����, �������������������Ѿ���ǰ��ϴ��(����˳��)
class CCombinationSampler :public CSampler
{
public:
	explicit CCombinationSampler(size_t numSamples);
	//explicit CCombinationSampler(const CEstimator& estimator);

	void Initialize(size_t numTotalSamples) override;
	size_t GetMaxNumSamples() const override;
	void Sample(std::vector<size_t>& sampledIndex) override;

private:
	std::vector<size_t> allIndex;
};











