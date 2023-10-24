#pragma once
#define EIGEN_USE_MKL_ALL
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <memory>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <any>
#include <Windows.h>

#include <boost/stacktrace.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>

#include "Base.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// 正数返回1, 负数返回-1, 0返回0
template <typename T>
int SignOfNumber(T val) noexcept
{
	return (T(0) < val) - (val < T(0));
}

// 将value限制在low和high的范围内
template <typename T>
T Clamp(const T& value, const T& low, const T& high) noexcept
{
	return std::max(low, std::min(value, high));
}

// 角度转弧度
inline double DegToRad(double deg) noexcept
{
	return deg * 0.0174532925199432954743716805978692718781530857086181640625;
}

// 弧度转角度
inline double RadToDeg(double rad) noexcept
{
	return rad * 57.29577951308232286464772187173366546630859375;
}

// 计算vector的中位数
template <typename T>
double GetMedian(const std::vector<T>& nums)
{
	CHECK(!nums.empty());

	const size_t midIndex = nums.size() / 2;
	std::vector<T> orderedNums = nums;
	std::nth_element(orderedNums.begin(), orderedNums.begin() + midIndex, orderedNums.end());

	if (nums.size() % 2 == 0)
	{
		const T midNum1 = orderedNums[midIndex];
		const T midNum2 = *std::max_element(orderedNums.begin(), orderedNums.begin() + midIndex);
		return 0.5 * midNum1 + 0.5 * midNum2;
	}
	return orderedNums[midIndex];
}

// 计算vector的平均数
template <typename T>
double GetMean(const std::vector<T>& nums)
{
	CHECK(!nums.empty());
	double sum = std::accumulate(nums.begin(), nums.end(), 0.0);
	return sum / nums.size();
}

// 计算vector的样本方差(无偏估计方差)
template <typename T>
double GetVariance(const std::vector<T>& nums)
{
	const double mean = GetMean(nums);
	double variance = 0;
	for (const auto num : nums)
	{
		const double diff = num - mean;
		variance += diff * diff;
	}
	return variance / (nums.size() - 1);
}

// 计算vector的标准差
template <typename T>
double GetStdDev(const std::vector<T>& nums)
{
	return std::sqrt(GetVariance(nums));
}

// 计算n个中选k个的组合数
size_t NChooseK(size_t n, size_t k);

// 生成n选k的组合. 注意: [first, last)中的元素必须按照std::less排序(升序)
template <typename Iterator>
bool NextCombination(Iterator first1, Iterator last1, Iterator first2, Iterator last2)
{
	if ((first1 == last1) || (first2 == last2))
	{
		return false;
	}
	Iterator m1 = last1, m2 = last2;
	m2--;
	while (--m1 != first1 && *m1 >= *m2) {}
	bool result = (m1 == first1) && (*first1 >= *m2);
	if (!result) 
	{
		while (first2 != m2 && *m1 >= *first2) 
		{
			first2++;
		}
		first1 = m1;
		std::iter_swap(first1, first2);
		first1++;
		first2++;
	}
	if ((first1 != last1) && (first2 != last2)) 
	{
		m1 = last1;
		m2 = first2;
		while ((m1 != first1) && (m2 != last2)) 
		{
			std::iter_swap(--m1, m2);
			m2++;
		}
		std::reverse(first1, m1);
		std::reverse(first1, last1);
		std::reverse(m2, last2);
		std::reverse(first2, last2);
	}
	return !result;
}

// 生成n选k的组合. 注意: [first, last)中的元素必须按照std::less排序(升序)
template <typename Iterator>
bool NextCombination(Iterator first, Iterator middle, Iterator last)
{
	return NextCombination(first, middle, middle, last);
}

// 将T1类型的值转换为T2类型，并在值超出T2类型的范围时进行截断
template <typename T1, typename T2>
T2 TruncateCast(T1 value)
{
	return static_cast<T2>(std::min(static_cast<T1>(std::numeric_limits<T2>::max()), std::max(static_cast<T1>(std::numeric_limits<T2>::min()), value)));
}

// 计算序列的第n百分位数
template <typename T>
T Percentile(const std::vector<T>& nums, double n)
{
	CHECK(!nums.empty());
	CHECK(n >= 0 && n <= 100);

	const int index = static_cast<int>(std::round(n / 100 * (nums.size() - 1)));
	const size_t percentileIndex = std::max(0, std::min(static_cast<int>(nums.size() - 1), index));

	std::vector<T> orderedNums = nums;
	std::nth_element(orderedNums.begin(), orderedNums.begin() + percentileIndex, orderedNums.end());

	return orderedNums[percentileIndex];
}

// 执行矩阵的RQ分解. 将矩阵A分解成上三角矩阵R和正交矩阵Q: A = R * Q
template <typename MatrixType>
void DecomposeMatrixRQ(const MatrixType& A, MatrixType& R, MatrixType& Q)
{
	const MatrixType A_flipud_transpose = A.transpose().rowwise().reverse().eval();

	const Eigen::HouseholderQR<MatrixType> QR(A_flipud_transpose);
	const MatrixType& Q0 = QR.householderQ();
	const MatrixType& R0 = QR.matrixQR();

	R = R0.transpose().colwise().reverse().eval();
	R = R.rowwise().reverse().eval();
	for (int i = 0; i < R.rows(); i++) 
	{
		for (int j = 0; j < R.cols() && (R.cols() - j) >(R.rows() - i); j++) 
		{
			R(i, j) = 0;
		}
	}

	Q = Q0.transpose().colwise().reverse().eval();

	if (Q.determinant() < 0) 
	{
		Q.row(1) *= -1.0;
		R.col(1) *= -1.0;
	}
}

// 移除多项式系数的前导零
Eigen::VectorXd RemoveLeadingZeros(const Eigen::VectorXd& coeffs);

// 移除多项式系数的尾随零
Eigen::VectorXd RemoveTrailingZeros(const Eigen::VectorXd& coeffs);

// 使用Horner方案计算多项式在x处的值, 多项式的各系数为coeffs
template <typename T>
T EvaluatePolynomial(const Eigen::VectorXd& coeffs, const T& x)
{
	T value = 0.0;
	for (Eigen::VectorXd::Index i = 0; i < coeffs.size(); i++) 
	{
		value = value * x + coeffs(i);
	}
	return value;
}

// 找到形如a*x+b=0的一次多项式的根, 将实部和虚部存储在real和imag中
bool FindLinearPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// 找到形如a*x^2+b*x+c=0的二次多项式的根, 将实部和虚部存储在real和imag中
bool FindQuadraticPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// 使用Durand-Kerner方法找到多项式的根. 该方法相对较快而且适用于任何多项式, 但可能在某些情况下不稳定或不准确
bool FindPolynomialRootsDurandKerner(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// 使用伴随矩阵方法找到多项式的根. 相比Durand-Kerner方法更慢, 但是更稳定和准确
// 参考文献: R. A. Horn & C. R. Johnson, Matrix Analysis. Cambridge, UK: Cambridge University Press, 1999, pp. 146-7.
bool FindPolynomialRootsCompanionMatrix(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag);

extern thread_local std::unique_ptr<std::mt19937> PRNG; // Mersenne Twister 19937 伪随机数生成器指针, 每个线程都有其自己的PRNG实例
extern int DefaultPRNGSeed; // 默认的PRNG种子

// 初始化PRNG. 如果种子seed为-1, 那么就会使用当前时间作为种子
void SetPRNGSeed(unsigned int seed = DefaultPRNGSeed);

// 生成在[min, max]范围的均匀分布的随机整数. 是无偏的, 线程安全的
template <typename T>
T GetRandomUniformInteger(T min, T max)
{
	if (!PRNG) 
	{
		SetPRNGSeed();
	}
	std::uniform_int_distribution<T> distribution(min, max);
	return distribution(*PRNG);
}

// 生成在[min, max]范围的均匀分布的随机实数. 是无偏的, 线程安全的
template <typename T>
T GetRandomUniformReal(T min, T max)
{
	if (!PRNG)
	{
		SetPRNGSeed();
	}
	std::uniform_real_distribution<T> distribution(min, max);
	return distribution(*PRNG);
}

// 生成一个满足平均值为mean, 标准差为stddev的高斯分布的随机实数. 是无偏的, 线程安全的
template <typename T>
T GetRandomGaussian(T mean, T stddev)
{
	if (PRNG == nullptr) 
	{
		SetPRNGSeed();
	}
	std::normal_distribution<T> distribution(mean, stddev);
	return distribution(*PRNG);
}

// 使用Fisher-Yates洗牌算法(也称Knuth洗牌算法)来随机排列一个向量的元素. numToShuffle表示要洗牌的元素数量, 这个参数允许仅洗牌向量中的一部分元素
template <typename T>
void Shuffle(size_t numToShuffle, std::vector<T>& nums)
{
	CHECK(numToShuffle <= nums.size());
	const size_t lastIndex = static_cast<size_t>(nums.size() - 1);
	for (size_t i = 0; i < numToShuffle; i++)
	{
		const size_t j = GetRandomUniformInteger<size_t>(i, lastIndex);
		std::swap(nums[i], nums[j]);
	}
}
