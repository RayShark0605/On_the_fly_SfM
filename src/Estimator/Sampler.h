#pragma once
#include <cstddef>
#include <vector>
#include <limits>

#include "../Base/Base.h"
#include "../Base/Math.h"
//#include "Estimator.h"

// 采样器. 这是一个抽象类
class CSampler
{
public:
	CSampler() = default;
	explicit CSampler(size_t numSamples);
	virtual ~CSampler() = default;

	// 初始化采样器. 在调用Sample之前必须调用这个方法.
	virtual void Initialize(size_t numTotalSamples) = 0;

	// 可以生成唯一样本的最大数量
	virtual size_t GetMaxNumSamples() const = 0;

	// 所有样本中采样numSamples个元素, 并将采样到的元素的索引存储在sampledIndex中
	virtual void Sample(std::vector<size_t>& sampledIndex) = 0;

	// 从一系列X中采样得到X_rand
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

	// 从一系列X中采样得到X_rand, 从一系列Y中采样得到Y_rand
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

// 随机采样器
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

// PROSAC方法(Progressive Sample Consensus)采样器. 参考文献: "Matching with PROSAC - Progressive Sample Consensus". Ondrej Chum and Matas, CVPR 2005.
// 每个线程都应该实例化一个单独的采样器, 而且质量更高的数据应该靠近列表的前端
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
	double T_n, T_n_p; // 在PROSAC算法中定义的与方程3相关的变量
};

// 用于基于RANSAC方法的随机采样器, 该采样器生成唯一的样本
// 每个线程应该实例化一个单独的采样器, 并且它假设输入数据已经提前被洗牌(打乱顺序)
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











