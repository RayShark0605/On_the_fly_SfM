#pragma once
#include <limits>
#include <vector>

// 模型的支持度
struct CSupport
{
	size_t numInliers = 0;                                    // 内点数量
	double residualSum = std::numeric_limits<double>::max();  // 内点残差总和(适用于CInlierSupportMeasurer)
	double score = std::numeric_limits<double>::max();        // MSAC得分, 定义为残差的截断和(适用于CMEstimatorSupportMeasurer)
};

// 模型支持度衡量器, 抽象类
class CSupportMeasurer
{
public:
	virtual CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept = 0; // 计算残差的支持度
	virtual bool Compare(const CSupport& support1, const CSupport& support2) const noexcept = 0;            // 比较support1是否比support2更好
};

// 通过计算内点inliers的数量和所有内点残差的总和来衡量模型的支持度. 如果一个支持度有更多的内点和更小的残差总和, 那么这个支持度就更好
class CInlierSupportMeasurer final :public CSupportMeasurer
{
public:
	CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept override;
	bool Compare(const CSupport& support1, const CSupport& support2) const noexcept override;
};

// 使用MSAC(M-estimator Sample and Consensus)方法中的数据适应度来衡量模型的支持度. 如果一个支持度有更小的MSAC得分, 那么这个支持度就更好
class CMEstimatorSupportMeasurer final : public CSupportMeasurer
{
public:
	CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept override;
	bool Compare(const CSupport& support1, const CSupport& support2) const noexcept override;
};

