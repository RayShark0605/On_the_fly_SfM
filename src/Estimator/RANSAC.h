#pragma once
#include <vector>

#include "../Base/Options.h"
#include "Sampler.h"
#include "SupportMeasurer.h"
#include "Estimator.h"

// RANSAC输出报告
template <typename ModelType>
struct CRANSACReport
{
	bool isSuccess = false;          // RANSAC估计是否成功
	size_t numTrials = 0;            // RANSAC迭代的次数
	CSupport support;                // 估计模型的支持度
	std::vector<char> inlierMask;    // 掩码, 是内点的样本就为true
	ModelType model;                 // RANSAC估计出的模型
};

// RANSAC算法, 用于从包含噪声和异常值的数据中稳健地估计一个模型
class CRANSAC
{
public:
	CRANSAC() = delete;
	CRANSAC(const CRANSACOptions& options, CEstimator* estimator, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr);
	virtual ~CRANSAC();

	// 确定在给定内点比例的情况下, 需要进行多少次尝试来保证以指定的置信度至少随机抽样出一个没有异常值的样本集
	// 
	// @param minNumSamples  估计模型所需的理论最少样本数
	// @param numInliers  内点的数量
	// @param numSamples  总样本数量
	// @param confidence  置信度(至少有一个样本是没有异常值的)
	// @param numTrialsMultiplier  尝试次数的乘法因子
	// 
	// @return  需要尝试的次数
	static size_t GetNumTrials(size_t minNumSamples, size_t numInliers, size_t numSamples, double confidence, double numTrialsMultiplier);
	size_t GetNumTrials(size_t numInliers, size_t numSamples, double confidence, double numTrialsMultiplier);

	// 使用RANSAC稳健地估计模型
	// 
	// @param X  自变量
	// @param Y  因变量
	// 
	// @return  RANSAC输出报告
	template <typename XType, typename YType, typename ModelType>
	CRANSACReport<ModelType> Estimate(const std::vector<XType>& X, const std::vector<YType>& Y)
	{
		// Step 1. 初始化. 检查输入向量大小, 初始化结果, 检查特殊情况
		Check(X.size() == Y.size());
		Check(estimator && supportMeasurer && sampler);

		const size_t numSamples = X.size();
		CRANSACReport<ModelType> report;
		if (numSamples < estimator->minNumSamples) // 如果当前样本数太少, 不足以用来估计模型, 则直接返回
		{
			return report;
		}

		// Step 2. RANSAC主循环
		bool isAbort = false;
		const double maxResidual = options.maxError * options.maxError; // 允许的最大残差(只有当残差不超过maxResidual, 才会被标记为内点)

		std::vector<double> residuals(numSamples); // 模型与所有数据点的残差
		std::vector<XType> X_rand(estimator->minNumSamples);
		std::vector<YType> Y_rand(estimator->minNumSamples);
		sampler->Initialize(numSamples); // 初始化采样器

		CSupport bestSupport; // 当前得到的最好的支持度
		ModelType bestModel;  // 当前得到的最好模型

		size_t maxNumTrials = std::min(options.maxNumTrials, sampler->GetMaxNumSamples()); // 确定最大迭代次数
		size_t dynamicMaxNumTrials = maxNumTrials;
		for (report.numTrials = 0; report.numTrials < maxNumTrials; report.numTrials++)
		{
			if (isAbort)
			{
				report.numTrials++;
				break;
			}

			sampler->GetSampleXY(X, Y, X_rand, Y_rand); // 从X和Y中以采样器内部的采样规则采样
			const std::vector<ModelType> sampledModels = AnyVec2TypeVec<ModelType>(estimator->Estimate(TypeVec2AnyVec(X_rand), TypeVec2AnyVec(Y_rand))); // 使用随机样本来估计模型

			for (const ModelType& sampledModel : sampledModels)
			{
				estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), sampledModel, residuals); // 对于每一个估计出的模型, 计算其与所有数据点的残差
				Check(residuals.size() == numSamples);

				const CSupport support = supportMeasurer->Evaluate(residuals, maxResidual); // 评估这个模型的支持度
				if (supportMeasurer->Compare(support, bestSupport)) // 如果新的支持度比当前最好的支持度更好, 就更新最好的支持度和最好的模型
				{
					bestSupport = support;
					bestModel = sampledModel;

					// 根据最好的支持度动态地更新最大迭代次数
					dynamicMaxNumTrials = GetNumTrials(estimator->minNumSamples, bestSupport.numInliers, numSamples, options.confidence, options.maxIterNumTrialsMultiplier);
				}

				// 如果达到动态更新后的最大迭代次数, 并且迭代次数超过了最小迭代次数，就提前终止算法
				if (report.numTrials >= dynamicMaxNumTrials && report.numTrials >= options.minNumTrials)
				{
					isAbort = true;
					break;
				}
			}
		}

		// Step 3. 结果收集与返回
		report.support = bestSupport;
		report.model = bestModel;
		if (report.support.numInliers < estimator->minNumSamples) // 如果找到的最佳模型的内点数少于最小样本数, 那么说明失败了
		{
			return report;
		}
		estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), report.model, residuals);
		Check(residuals.size() == numSamples);

		report.inlierMask.resize(numSamples);
		for (size_t i = 0; i < residuals.size(); i++) // 判断每个样本是否为内点
		{
			report.inlierMask[i] = (residuals[i] <= maxResidual);
		}
		report.isSuccess = true;
		return report;
	}

protected:
	CRANSACOptions options;
	CEstimator* estimator;
	CSupportMeasurer* supportMeasurer;
	CSampler* sampler;

	bool isSupportMeasurerInitialNull = false;
	bool isSamplerInitialNull = false;
};


// LORANSAC(增加了局部优化的RANSAC)实现. 参考文献: "Locally Optimized RANSAC" Ondrej Chum, Jiri Matas, Josef Kittler, DAGM 2003.
// 在迭代过程中, 一旦找到一个更好的模型, 尝试通过局部搜索, 增加更多内点来进一步优化它
class CLORANSAC: public CRANSAC
{
public:
	CLORANSAC() = delete;
	CLORANSAC(const CRANSACOptions& options, CEstimator* estimator, CEstimator* localEstimator, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr);

	// 使用LORANSAC稳健地估计模型
	// 
	// @param X  自变量
	// @param Y  因变量
	// 
	// @return  RANSAC输出报告
	template <typename XType, typename YType, typename ModelType>
	CRANSACReport<ModelType> Estimate(const std::vector<XType>& X, const std::vector<YType>& Y)
	{
		// Step 1. 初始化. 检查输入向量大小, 初始化结果, 检查特殊情况
		Check(X.size() == Y.size());
		Check(estimator && localEstimator && supportMeasurer && sampler);

		const size_t numSamples = X.size();
		CRANSACReport<ModelType> report;
		if (numSamples < estimator->minNumSamples) // 如果当前样本数太少, 不足以用来估计模型, 则直接返回
		{
			return report;
		}

		// Step 2. RANSAC主循环
		bool isAbort = false;
		bool bestModelIsLocal = false;
		const double maxResidual = options.maxError * options.maxError; // 允许的最大残差(只有当残差不超过maxResidual, 才会被标记为内点)

		std::vector<double> residuals; // 模型与所有数据点的残差
		std::vector<double> bestLocalResiduals; // 模型与所有数据点的残差
		std::vector<XType> XInlier;
		std::vector<YType> YInlier;
		std::vector<XType> X_rand(estimator->minNumSamples);
		std::vector<YType> Y_rand(estimator->minNumSamples);
		sampler->Initialize(numSamples); // 初始化采样器

		CSupport bestSupport; // 当前得到的最好的支持度
		ModelType bestModel;  // 当前得到的最好模型

		size_t maxNumTrials = std::min(options.maxNumTrials, sampler->GetMaxNumSamples()); // 确定最大迭代次数
		size_t dynamicMaxNumTrials = maxNumTrials;
		for (report.numTrials = 0; report.numTrials < maxNumTrials; report.numTrials++)
		{
			if (isAbort)
			{
				report.numTrials++;
				break;
			}

			sampler->GetSampleXY(X, Y, X_rand, Y_rand); // 从X和Y中以采样器内部的采样规则采样
			const std::vector<ModelType> sampledModels = AnyVec2TypeVec<ModelType>(estimator->Estimate(TypeVec2AnyVec(std::move(X_rand)), TypeVec2AnyVec(std::move(Y_rand)))); // 使用随机样本来估计模型
			for (const ModelType& sampledModel : sampledModels)
			{
				estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), sampledModel, residuals); // 对于每一个估计出的模型, 计算其与所有数据点的残差
				Check(residuals.size() == numSamples);

				const CSupport support = supportMeasurer->Evaluate(residuals, maxResidual); // 评估这个模型的支持度

				if (supportMeasurer->Compare(support, bestSupport)) // 如果新的支持度比当前最好的支持度更好, 就做局部优化
				{
					bestSupport = support;
					bestModel = sampledModel;
					bestModelIsLocal = false;

					// 根据内点来局部估计更优模型
					if (support.numInliers > estimator->minNumSamples && support.numInliers >= localEstimator->minNumSamples)
					{
						// 迭代式局部优化来扩大内点集
						const size_t maxLocalTrials = 10;
						for (size_t localNumTrials = 0; localNumTrials < maxLocalTrials; localNumTrials++)
						{
							XInlier.clear();
							YInlier.clear();
							XInlier.reserve(numSamples);
							YInlier.reserve(numSamples);
							for (size_t i = 0; i < residuals.size(); i++)
							{
								if (residuals[i] <= maxResidual)
								{
									XInlier.push_back(X[i]);
									YInlier.push_back(Y[i]);
								}
							}
							const std::vector<ModelType> localModels = AnyVec2TypeVec<ModelType>(localEstimator->Estimate(TypeVec2AnyVec(std::move(XInlier)), TypeVec2AnyVec(std::move(YInlier))));
							const size_t preBestNumInliers = bestSupport.numInliers;
							for (const ModelType& localModel : localModels)
							{
								localEstimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), localModel, residuals);
								Check(residuals.size() == numSamples);
								const CSupport localSupport = supportMeasurer->Evaluate(residuals, maxResidual);

								// 检查局部优化模型是否更优
								if (supportMeasurer->Compare(localSupport, bestSupport))
								{
									bestSupport = localSupport;
									bestModel = localModel;
									bestModelIsLocal = true;
									std::swap(residuals, bestLocalResiduals); // 交换残差
								}
							}
							// 只有当内点集变多了从而有机会进一步优化时, 才继续迭代
							if (bestSupport.numInliers <= preBestNumInliers)
							{
								break;
							}

							//把残差再交换回来, 这样就可以在下一次局部优化的迭代中提取出最佳的内点集
							std::swap(residuals, bestLocalResiduals);
						}
					}
					dynamicMaxNumTrials = GetNumTrials(estimator->minNumSamples, bestSupport.numInliers, numSamples, options.confidence, options.maxIterNumTrialsMultiplier);
				}
				if (report.numTrials >= dynamicMaxNumTrials && report.numTrials >= options.minNumTrials)
				{
					isAbort = true;
					break;
				}
			}
		}

		// Step 3. 结果收集与返回
		report.support = bestSupport;
		report.model = bestModel;
		if (report.support.numInliers < estimator->minNumSamples) // 如果找到的最佳模型的内点数少于最小样本数, 那么说明失败了
		{
			return report;
		}

		// 这将对最佳模型的残差进行两次计算, 但避免了对每个评估模型都复制和填充内点掩码, 这种方法其实更快
		if (bestModelIsLocal)
		{
			localEstimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), report.model, residuals);
		}
		else
		{
			estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), report.model, residuals);
		}
		Check(residuals.size() == numSamples);
		report.inlierMask.resize(numSamples);

		for (size_t i = 0; i < residuals.size(); i++) // 判断每个样本是否为内点
		{
			report.inlierMask[i] = (residuals[i] <= maxResidual);
		}
		report.isSuccess = true;
		return report;
	}

private:
	CEstimator* localEstimator;
};
