#include "Math.h"
#include "Base.h"

thread_local std::unique_ptr<std::mt19937> PRNG;
int DefaultPRNGSeed = 0;

size_t NChooseK(size_t n, size_t k)
{
	if (k == 0 || k == n) return 1;
	if (k > n) return 0;
	if (k > n - k) k = n - k;

	size_t result = 1;
	for (size_t i = 0; i < k; i++)
	{
		result *= (n - i);
		result /= (i + 1);
	}
	return result;
}
Eigen::VectorXd RemoveLeadingZeros(const Eigen::VectorXd& coeffs)
{
	Eigen::VectorXd::Index numZeros = 0;
	for (; numZeros < coeffs.size(); numZeros++)
	{
		if (coeffs(numZeros) != 0)
		{
			break;
		}
	}
	return coeffs.tail(coeffs.size() - numZeros);
}
Eigen::VectorXd RemoveTrailingZeros(const Eigen::VectorXd& coeffs)
{
	Eigen::VectorXd::Index numZeros = 0;
	for (; numZeros < coeffs.size(); numZeros++)
	{
		if (coeffs(coeffs.size() - 1 - numZeros) != 0)
		{
			break;
		}
	}
	return coeffs.head(coeffs.size() - numZeros);
}
bool FindLinearPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag)
{
	CHECK(coeffs.size() == 2);
	if (coeffs(0) == 0)
	{
		return false;
	}
	real.resize(1);
	real(0) = -coeffs(1) / coeffs(0);

	imag.resize(1);
	imag(0) = 0;
	return true;
}
bool FindQuadraticPolynomialRoots(const Eigen::VectorXd& coeffs, Eigen::VectorXd& real, Eigen::VectorXd& imag)
{
	CHECK(coeffs.size() == 3);
	const double a = coeffs(0);
	if (a == 0)
	{
		return FindLinearPolynomialRoots(coeffs.tail(2), real, imag);
	}
	const double b = coeffs(1);
	const double c = coeffs(2);
	if (b == 0 && c == 0)
	{
		real.resize(1);
		real(0) = 0;

		imag.resize(1);
		imag(0) = 0;

		return true;
	}
	const double d = b * b - 4 * a * c;
	if (d >= 0)
	{
		const double sqrt_d = std::sqrt(d);
		real.resize(2);
		if (b >= 0)
		{
			real(0) = (-b - sqrt_d) / (2 * a);
			real(1) = (2 * c) / (-b - sqrt_d);
		}
		else
		{
			real(0) = (2 * c) / (-b + sqrt_d);
			real(1) = (-b + sqrt_d) / (2 * a);
		}
		imag.resize(2);
		imag.setZero();
	}
	else
	{
		real.resize(2);
		real.setConstant(-b / (2 * a));

		imag.resize(2);
		imag(0) = std::sqrt(-d) / (2 * a);
		imag(1) = -imag(0);
	}
	return true;
}
bool FindPolynomialRootsDurandKerner(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag)
{
	CHECK(coeffsAll.size() >= 2);

	// Step 1. 输入检查和预处理, 去除前导零, 确定多项式的度
	const Eigen::VectorXd coeffs = RemoveLeadingZeros(coeffsAll);
	const int degree = coeffs.size() - 1;

	// Step 2. 特殊情况处理
	if (degree <= 0)
	{
		return false;
	}
	if (degree == 1)
	{
		return FindLinearPolynomialRoots(coeffs, real, imag);
	}
	if (degree == 2)
	{
		return FindQuadraticPolynomialRoots(coeffs, real, imag);
	}

	// Step 3. 初始化根
	Eigen::VectorXcd roots(degree);
	roots(degree - 1) = std::complex<double>(1, 0);
	for (int i = degree - 2; i >= 0; i--)
	{
		roots(i) = roots(i + 1) * std::complex<double>(1, 1);
	}

	// Step 4. 迭代求解
	const int maxNumIterations = 100;
	const double maxRootChange = 1e-10;
	for (int iter = 0; iter < maxNumIterations; iter++)
	{
		double curMaxRootChange = 0.0;
		for (int i = 0; i < degree; i++)
		{
			const std::complex<double> root_i = roots(i);
			std::complex<double> numerator = coeffs[0];
			std::complex<double> denominator = coeffs[0];
			for (int j = 0; j < degree; j++)
			{
				numerator = numerator * root_i + coeffs[j + 1];
				if (i != j)
				{
					denominator = denominator * (root_i - roots(j));
				}
			}
			const std::complex<double> root_i_Change = numerator / denominator;
			roots(i) = root_i - root_i_Change;
			curMaxRootChange = std::max(curMaxRootChange, std::abs(root_i_Change.real()));
			curMaxRootChange = std::max(curMaxRootChange, std::abs(root_i_Change.imag()));
		}
		if (curMaxRootChange < maxRootChange)
		{
			break;
		}
	}

	// Step 5. 输出结果
	real.resize(degree);
	real = roots.real();

	imag.resize(degree);
	imag = roots.imag();

	return true;
}
bool FindPolynomialRootsCompanionMatrix(const Eigen::VectorXd& coeffsAll, Eigen::VectorXd& real, Eigen::VectorXd& imag)
{
	CHECK(coeffsAll.size() >= 2);

	// Step 1. 输入检查和预处理, 去除前导零, 确定多项式的度
	Eigen::VectorXd coeffs = RemoveLeadingZeros(coeffsAll);
	const int degree = coeffs.size() - 1;

	// Step 2. 特殊情况处理
	if (degree <= 0)
	{
		return false;
	}
	if (degree == 1)
	{
		return FindLinearPolynomialRoots(coeffs, real, imag);
	}
	if (degree == 2)
	{
		return FindQuadraticPolynomialRoots(coeffs, real, imag);
	}

	// Step 3. 移除尾随零并检查是否只有零是解
	coeffs = RemoveTrailingZeros(coeffs);
	if (coeffs.size() == 1)
	{
		real.resize(1);
		real(0) = 0;
		imag.resize(1);
		imag(0) = 0;
		return true;
	}

	// Step 4. 填充伴随矩阵C
	Eigen::MatrixXd C(coeffs.size() - 1, coeffs.size() - 1);
	C.setZero();
	for (Eigen::MatrixXd::Index i = 1; i < C.rows(); i++)
	{
		C(i, i - 1) = 1;
	}
	C.row(0) = -coeffs.tail(coeffs.size() - 1) / coeffs(0);

	// Step 5. 使用EigenSolver来求解伴随矩阵的特征值
	Eigen::EigenSolver<Eigen::MatrixXd> solver(C, false);
	if (solver.info() != Eigen::Success)
	{
		return false;
	}

	// Step 6. 处理尾随零和输出结果. 如果存在尾随零, 意味着零也是一个解
	const int effectiveDegree = coeffs.size() - 1 < degree ? coeffs.size() : coeffs.size() - 1;

	real.resize(effectiveDegree);
	real.head(coeffs.size() - 1) = solver.eigenvalues().real();
	if (effectiveDegree > coeffs.size() - 1)
	{
		real(real.size() - 1) = 0;
	}

	imag.resize(effectiveDegree);
	imag.head(coeffs.size() - 1) = solver.eigenvalues().imag();
	if (effectiveDegree > coeffs.size() - 1)
	{
		imag(imag.size() - 1) = 0;
	}

	return true;
}
void SetPRNGSeed(unsigned int seed)
{
	PRNG = std::make_unique<std::mt19937>(seed);
	static std::mutex mutex;
	std::unique_lock<std::mutex> lock(mutex);
	srand(seed);
}




