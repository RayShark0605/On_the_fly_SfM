#include <numeric>
#include <limits>

#include "Sampler.h"

using namespace std;

CSampler::CSampler(size_t numSamples) :numSamples(numSamples)
{

}

CRandomSampler::CRandomSampler(size_t numSamples) :CSampler(numSamples)
{
	numTotalSamples = 0;
}
//CRandomSampler::CRandomSampler(const CEstimator& estimator): CSampler(estimator.minNumSamples)
//{
//	numTotalSamples = 0;
//}
void CRandomSampler::Initialize(size_t numTotalSamples)
{
	Check(numSamples <= numTotalSamples);
	this->numTotalSamples = numTotalSamples;
	allIndex.resize(numTotalSamples);
	iota(allIndex.begin(), allIndex.end(), 0);
}
size_t CRandomSampler::GetMaxNumSamples() const
{
	return numeric_limits<size_t>::max();
}
void CRandomSampler::Sample(vector<size_t>& sampledIndex)
{
	Shuffle(numSamples, allIndex);
	sampledIndex.resize(numSamples);
	for (size_t i = 0; i < numSamples; i++)
	{
		sampledIndex[i] = allIndex[i];
	}
}

CProgressiveSampler::CProgressiveSampler(size_t numSamples) :CSampler(numSamples)
{
	numTotalSamples = 0;
	t = 0;
	n = 0;
	T_n = 0;
	T_n_p = 0;
}
//CProgressiveSampler::CProgressiveSampler(const CEstimator& estimator) :CSampler(estimator.minNumSamples)
//{
//	numTotalSamples = 0;
//	t = 0;
//	n = 0;
//	T_n = 0;
//	T_n_p = 0;
//}
void CProgressiveSampler::Initialize(size_t numTotalSamples)
{
	Check(numSamples <= numTotalSamples);
	this->numTotalSamples = numTotalSamples;
	t = 0;
	n = numSamples;

	T_n = 200000; // 在PROSAC表现得像RANSAC之前的迭代次数, 论文推荐值
	T_n_p = 1;
	for (size_t i = 0; i < numSamples; i++)
	{
		T_n *= static_cast<double>(numSamples - i) / (numTotalSamples - i);
	}
}
size_t CProgressiveSampler::GetMaxNumSamples() const
{
	return numeric_limits<size_t>::max();
}
void CProgressiveSampler::Sample(vector<size_t>& sampledIndex)
{
	t++;
	sampledIndex.clear();
	sampledIndex.reserve(numSamples);
	if (abs(t - T_n_p) < 1e-5 && n < numTotalSamples)
	{
		const double T_n_plus_1 = T_n * (n + 1.0) / (n + 1.0 - numSamples);
		T_n_p += ceil(T_n_plus_1 - T_n);
		T_n = T_n_plus_1;
		n++;
	}

	Check(n > 1 && numSamples > 0);
	size_t numRandomSamples = numSamples;
	size_t maxRandomSampleIndex = n - 1;
	if (T_n_p >= t) 
	{
		numRandomSamples--;
		maxRandomSampleIndex--;
	}

	unordered_set<size_t> set;
	set.reserve(numSamples);
	for (size_t i = 0; i < numRandomSamples; i++)
	{
		while (true) 
		{
			const size_t randomIndex = GetRandomUniformInteger<size_t>(0, maxRandomSampleIndex);
			if (set.find(randomIndex) == set.end())
			{
				sampledIndex.push_back(randomIndex);
				set.insert(randomIndex);
				break;
			}
		}
	}
	if (T_n_p >= t) 
	{
		sampledIndex.push_back(n);
	}
}

CCombinationSampler::CCombinationSampler(size_t numSamples) :CSampler(numSamples)
{
	numTotalSamples = 0;
	allIndex.clear();
}
//CCombinationSampler::CCombinationSampler(const CEstimator& estimator) :CSampler(estimator.minNumSamples)
//{
//	numTotalSamples = 0;
//	allIndex.clear();
//}
void CCombinationSampler::Initialize(size_t numTotalSamples)
{
	Check(numSamples <= numTotalSamples);
	this->numTotalSamples = numTotalSamples;
	allIndex.resize(numTotalSamples);
	iota(allIndex.begin(), allIndex.end(), 0);
}
size_t CCombinationSampler::GetMaxNumSamples() const
{
	return NChooseK(allIndex.size(), numSamples);
}
void CCombinationSampler::Sample(vector<size_t>& sampledIndex)
{
	sampledIndex.resize(numSamples);
	for (size_t i = 0; i < numSamples; i++)
	{
		sampledIndex[i] = allIndex[i];
	}
	if (!NextCombination(allIndex.begin(), allIndex.begin() + numSamples, allIndex.end()))
	{
		iota(allIndex.begin(), allIndex.end(), 0); // 达到了所有可能的组合, 因此重置为原始状态
	}
}


