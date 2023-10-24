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

// ��������1, ��������-1, 0����0
template <typename T>
int SignOfNumber(T val) noexcept
{
	return (T(0) < val) - (val < T(0));
}

// ��value������low��high�ķ�Χ��
template <typename T>
T Clamp(const T& value, const T& low, const T& high) noexcept
{
	return std::max(low, std::min(value, high));
}

// �Ƕ�ת����
inline double DegToRad(double deg) noexcept
{
	return deg * 0.0174532925199432954743716805978692718781530857086181640625;
}

// ����ת�Ƕ�
inline double RadToDeg(double rad) noexcept
{
	return rad * 57.29577951308232286464772187173366546630859375;
}

// ����vector����λ��
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

// ����vector��ƽ����
template <typename T>
double GetMean(const std::vector<T>& nums)
{
	CHECK(!nums.empty());
	double sum = std::accumulate(nums.begin(), nums.end(), 0.0);
	return sum / nums.size();
}

// ����vector����������(��ƫ���Ʒ���)
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

// ����vector�ı�׼��
template <typename T>
double GetStdDev(const std::vector<T>& nums)
{
	return std::sqrt(GetVariance(nums));
}

// ����n����ѡk���������
size_t NChooseK(size_t n, size_t k);

// ����nѡk�����. ע��: [first, last)�е�Ԫ�ر��밴��std::less����(����)
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

// ����nѡk�����. ע��: [first, last)�е�Ԫ�ر��밴��std::less����(����)
template <typename Iterator>
bool NextCombination(Iterator first, Iterator middle, Iterator last)
{
	return NextCombination(first, middle, middle, last);
}

// ��T1���͵�ֵת��ΪT2���ͣ�����ֵ����T2���͵ķ�Χʱ���нض�
template <typename T1, typename T2>
T2 TruncateCast(T1 value)
{
	return static_cast<T2>(std::min(static_cast<T1>(std::numeric_limits<T2>::max()), std::max(static_cast<T1>(std::numeric_limits<T2>::min()), value)));
}

// �������еĵ�n�ٷ�λ��
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

// ִ�о����RQ�ֽ�. ������A�ֽ�������Ǿ���R����������Q: A = R * Q
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

// �Ƴ�����ʽϵ����ǰ����
Eigen::VectorXd RemoveLeadingZeros(const Eigen::VectorXd& coeffs);

// �Ƴ�����ʽϵ����β����
Eigen::VectorXd RemoveTrailingZeros(const Eigen::VectorXd& coeffs);

// ʹ��Horner�����������ʽ��x����ֵ, ����ʽ�ĸ�ϵ��Ϊcoeffs
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

// �ҵ�����a*x+b=0��һ�ζ���ʽ�ĸ�, ��ʵ�����鲿�洢��real��imag��
bool FindLinearPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// �ҵ�����a*x^2+b*x+c=0�Ķ��ζ���ʽ�ĸ�, ��ʵ�����鲿�洢��real��imag��
bool FindQuadraticPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// ʹ��Durand-Kerner�����ҵ�����ʽ�ĸ�. �÷�����ԽϿ�����������κζ���ʽ, ��������ĳЩ����²��ȶ���׼ȷ
bool FindPolynomialRootsDurandKerner(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag);

// ʹ�ð�����󷽷��ҵ�����ʽ�ĸ�. ���Durand-Kerner��������, ���Ǹ��ȶ���׼ȷ
// �ο�����: R. A. Horn & C. R. Johnson, Matrix Analysis. Cambridge, UK: Cambridge University Press, 1999, pp. 146-7.
bool FindPolynomialRootsCompanionMatrix(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag);

extern thread_local std::unique_ptr<std::mt19937> PRNG; // Mersenne Twister 19937 α�����������ָ��, ÿ���̶߳������Լ���PRNGʵ��
extern int DefaultPRNGSeed; // Ĭ�ϵ�PRNG����

// ��ʼ��PRNG. �������seedΪ-1, ��ô�ͻ�ʹ�õ�ǰʱ����Ϊ����
void SetPRNGSeed(unsigned int seed = DefaultPRNGSeed);

// ������[min, max]��Χ�ľ��ȷֲ����������. ����ƫ��, �̰߳�ȫ��
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

// ������[min, max]��Χ�ľ��ȷֲ������ʵ��. ����ƫ��, �̰߳�ȫ��
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

// ����һ������ƽ��ֵΪmean, ��׼��Ϊstddev�ĸ�˹�ֲ������ʵ��. ����ƫ��, �̰߳�ȫ��
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

// ʹ��Fisher-Yatesϴ���㷨(Ҳ��Knuthϴ���㷨)���������һ��������Ԫ��. numToShuffle��ʾҪϴ�Ƶ�Ԫ������, ������������ϴ�������е�һ����Ԫ��
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
