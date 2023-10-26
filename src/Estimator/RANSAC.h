#pragma once
#include <vector>

#include "../Base/Options.h"
#include "Sampler.h"
#include "SupportMeasurer.h"
#include "Estimator.h"

// RANSAC�������
template <typename ModelType>
struct CRANSACReport
{
	bool isSuccess = false;          // RANSAC�����Ƿ�ɹ�
	size_t numTrials = 0;            // RANSAC�����Ĵ���
	CSupport support;                // ����ģ�͵�֧�ֶ�
	std::vector<char> inlierMask;    // ����, ���ڵ��������Ϊtrue
	ModelType model;                 // RANSAC���Ƴ���ģ��
};

// RANSAC�㷨, ���ڴӰ����������쳣ֵ���������Ƚ��ع���һ��ģ��
class CRANSAC
{
public:
	CRANSAC() = delete;
	CRANSAC(const CRANSACOptions& options, CEstimator* estimator, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr);
	virtual ~CRANSAC();

	// ȷ���ڸ����ڵ�����������, ��Ҫ���ж��ٴγ�������֤��ָ�������Ŷ��������������һ��û���쳣ֵ��������
	// 
	// @param minNumSamples  ����ģ���������������������
	// @param numInliers  �ڵ������
	// @param numSamples  ����������
	// @param confidence  ���Ŷ�(������һ��������û���쳣ֵ��)
	// @param numTrialsMultiplier  ���Դ����ĳ˷�����
	// 
	// @return  ��Ҫ���ԵĴ���
	static size_t GetNumTrials(size_t minNumSamples, size_t numInliers, size_t numSamples, double confidence, double numTrialsMultiplier);
	size_t GetNumTrials(size_t numInliers, size_t numSamples, double confidence, double numTrialsMultiplier);

	// ʹ��RANSAC�Ƚ��ع���ģ��
	// 
	// @param X  �Ա���
	// @param Y  �����
	// 
	// @return  RANSAC�������
	template <typename XType, typename YType, typename ModelType>
	CRANSACReport<ModelType> Estimate(const std::vector<XType>& X, const std::vector<YType>& Y)
	{
		// Step 1. ��ʼ��. �������������С, ��ʼ�����, ����������
		Check(X.size() == Y.size());
		Check(estimator && supportMeasurer && sampler);

		const size_t numSamples = X.size();
		CRANSACReport<ModelType> report;
		if (numSamples < estimator->minNumSamples) // �����ǰ������̫��, ��������������ģ��, ��ֱ�ӷ���
		{
			return report;
		}

		// Step 2. RANSAC��ѭ��
		bool isAbort = false;
		const double maxResidual = options.maxError * options.maxError; // ��������в�(ֻ�е��в����maxResidual, �Żᱻ���Ϊ�ڵ�)

		std::vector<double> residuals(numSamples); // ģ�����������ݵ�Ĳв�
		std::vector<XType> X_rand(estimator->minNumSamples);
		std::vector<YType> Y_rand(estimator->minNumSamples);
		sampler->Initialize(numSamples); // ��ʼ��������

		CSupport bestSupport; // ��ǰ�õ�����õ�֧�ֶ�
		ModelType bestModel;  // ��ǰ�õ������ģ��

		size_t maxNumTrials = std::min(options.maxNumTrials, sampler->GetMaxNumSamples()); // ȷ������������
		size_t dynamicMaxNumTrials = maxNumTrials;
		for (report.numTrials = 0; report.numTrials < maxNumTrials; report.numTrials++)
		{
			if (isAbort)
			{
				report.numTrials++;
				break;
			}

			sampler->GetSampleXY(X, Y, X_rand, Y_rand); // ��X��Y���Բ������ڲ��Ĳ����������
			const std::vector<ModelType> sampledModels = AnyVec2TypeVec<ModelType>(estimator->Estimate(TypeVec2AnyVec(X_rand), TypeVec2AnyVec(Y_rand))); // ʹ���������������ģ��

			for (const ModelType& sampledModel : sampledModels)
			{
				estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), sampledModel, residuals); // ����ÿһ�����Ƴ���ģ��, ���������������ݵ�Ĳв�
				Check(residuals.size() == numSamples);

				const CSupport support = supportMeasurer->Evaluate(residuals, maxResidual); // �������ģ�͵�֧�ֶ�
				if (supportMeasurer->Compare(support, bestSupport)) // ����µ�֧�ֶȱȵ�ǰ��õ�֧�ֶȸ���, �͸�����õ�֧�ֶȺ���õ�ģ��
				{
					bestSupport = support;
					bestModel = sampledModel;

					// ������õ�֧�ֶȶ�̬�ظ�������������
					dynamicMaxNumTrials = GetNumTrials(estimator->minNumSamples, bestSupport.numInliers, numSamples, options.confidence, options.maxIterNumTrialsMultiplier);
				}

				// ����ﵽ��̬���º������������, ���ҵ���������������С��������������ǰ��ֹ�㷨
				if (report.numTrials >= dynamicMaxNumTrials && report.numTrials >= options.minNumTrials)
				{
					isAbort = true;
					break;
				}
			}
		}

		// Step 3. ����ռ��뷵��
		report.support = bestSupport;
		report.model = bestModel;
		if (report.support.numInliers < estimator->minNumSamples) // ����ҵ������ģ�͵��ڵ���������С������, ��ô˵��ʧ����
		{
			return report;
		}
		estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), report.model, residuals);
		Check(residuals.size() == numSamples);

		report.inlierMask.resize(numSamples);
		for (size_t i = 0; i < residuals.size(); i++) // �ж�ÿ�������Ƿ�Ϊ�ڵ�
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


// LORANSAC(�����˾ֲ��Ż���RANSAC)ʵ��. �ο�����: "Locally Optimized RANSAC" Ondrej Chum, Jiri Matas, Josef Kittler, DAGM 2003.
// �ڵ���������, һ���ҵ�һ�����õ�ģ��, ����ͨ���ֲ�����, ���Ӹ����ڵ�����һ���Ż���
class CLORANSAC: public CRANSAC
{
public:
	CLORANSAC() = delete;
	CLORANSAC(const CRANSACOptions& options, CEstimator* estimator, CEstimator* localEstimator, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr);

	// ʹ��LORANSAC�Ƚ��ع���ģ��
	// 
	// @param X  �Ա���
	// @param Y  �����
	// 
	// @return  RANSAC�������
	template <typename XType, typename YType, typename ModelType>
	CRANSACReport<ModelType> Estimate(const std::vector<XType>& X, const std::vector<YType>& Y)
	{
		// Step 1. ��ʼ��. �������������С, ��ʼ�����, ����������
		Check(X.size() == Y.size());
		Check(estimator && localEstimator && supportMeasurer && sampler);

		const size_t numSamples = X.size();
		CRANSACReport<ModelType> report;
		if (numSamples < estimator->minNumSamples) // �����ǰ������̫��, ��������������ģ��, ��ֱ�ӷ���
		{
			return report;
		}

		// Step 2. RANSAC��ѭ��
		bool isAbort = false;
		bool bestModelIsLocal = false;
		const double maxResidual = options.maxError * options.maxError; // ��������в�(ֻ�е��в����maxResidual, �Żᱻ���Ϊ�ڵ�)

		std::vector<double> residuals; // ģ�����������ݵ�Ĳв�
		std::vector<double> bestLocalResiduals; // ģ�����������ݵ�Ĳв�
		std::vector<XType> XInlier;
		std::vector<YType> YInlier;
		std::vector<XType> X_rand(estimator->minNumSamples);
		std::vector<YType> Y_rand(estimator->minNumSamples);
		sampler->Initialize(numSamples); // ��ʼ��������

		CSupport bestSupport; // ��ǰ�õ�����õ�֧�ֶ�
		ModelType bestModel;  // ��ǰ�õ������ģ��

		size_t maxNumTrials = std::min(options.maxNumTrials, sampler->GetMaxNumSamples()); // ȷ������������
		size_t dynamicMaxNumTrials = maxNumTrials;
		for (report.numTrials = 0; report.numTrials < maxNumTrials; report.numTrials++)
		{
			if (isAbort)
			{
				report.numTrials++;
				break;
			}

			sampler->GetSampleXY(X, Y, X_rand, Y_rand); // ��X��Y���Բ������ڲ��Ĳ����������
			const std::vector<ModelType> sampledModels = AnyVec2TypeVec<ModelType>(estimator->Estimate(TypeVec2AnyVec(std::move(X_rand)), TypeVec2AnyVec(std::move(Y_rand)))); // ʹ���������������ģ��
			for (const ModelType& sampledModel : sampledModels)
			{
				estimator->Residuals(TypeVec2AnyVec(X), TypeVec2AnyVec(Y), sampledModel, residuals); // ����ÿһ�����Ƴ���ģ��, ���������������ݵ�Ĳв�
				Check(residuals.size() == numSamples);

				const CSupport support = supportMeasurer->Evaluate(residuals, maxResidual); // �������ģ�͵�֧�ֶ�

				if (supportMeasurer->Compare(support, bestSupport)) // ����µ�֧�ֶȱȵ�ǰ��õ�֧�ֶȸ���, �����ֲ��Ż�
				{
					bestSupport = support;
					bestModel = sampledModel;
					bestModelIsLocal = false;

					// �����ڵ����ֲ����Ƹ���ģ��
					if (support.numInliers > estimator->minNumSamples && support.numInliers >= localEstimator->minNumSamples)
					{
						// ����ʽ�ֲ��Ż��������ڵ㼯
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

								// ���ֲ��Ż�ģ���Ƿ����
								if (supportMeasurer->Compare(localSupport, bestSupport))
								{
									bestSupport = localSupport;
									bestModel = localModel;
									bestModelIsLocal = true;
									std::swap(residuals, bestLocalResiduals); // �����в�
								}
							}
							// ֻ�е��ڵ㼯����˴Ӷ��л����һ���Ż�ʱ, �ż�������
							if (bestSupport.numInliers <= preBestNumInliers)
							{
								break;
							}

							//�Ѳв��ٽ�������, �����Ϳ�������һ�ξֲ��Ż��ĵ�������ȡ����ѵ��ڵ㼯
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

		// Step 3. ����ռ��뷵��
		report.support = bestSupport;
		report.model = bestModel;
		if (report.support.numInliers < estimator->minNumSamples) // ����ҵ������ģ�͵��ڵ���������С������, ��ô˵��ʧ����
		{
			return report;
		}

		// �⽫�����ģ�͵Ĳв�������μ���, �������˶�ÿ������ģ�Ͷ����ƺ�����ڵ�����, ���ַ�����ʵ����
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

		for (size_t i = 0; i < residuals.size(); i++) // �ж�ÿ�������Ƿ�Ϊ�ڵ�
		{
			report.inlierMask[i] = (residuals[i] <= maxResidual);
		}
		report.isSuccess = true;
		return report;
	}

private:
	CEstimator* localEstimator;
};
