#pragma once
#include <limits>
#include <vector>

// ģ�͵�֧�ֶ�
struct CSupport
{
	size_t numInliers = 0;                                    // �ڵ�����
	double residualSum = std::numeric_limits<double>::max();  // �ڵ�в��ܺ�(������CInlierSupportMeasurer)
	double score = std::numeric_limits<double>::max();        // MSAC�÷�, ����Ϊ�в�ĽضϺ�(������CMEstimatorSupportMeasurer)
};

// ģ��֧�ֶȺ�����, ������
class CSupportMeasurer
{
public:
	virtual CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept = 0; // ����в��֧�ֶ�
	virtual bool Compare(const CSupport& support1, const CSupport& support2) const noexcept = 0;            // �Ƚ�support1�Ƿ��support2����
};

// ͨ�������ڵ�inliers�������������ڵ�в���ܺ�������ģ�͵�֧�ֶ�. ���һ��֧�ֶ��и�����ڵ�͸�С�Ĳв��ܺ�, ��ô���֧�ֶȾ͸���
class CInlierSupportMeasurer final :public CSupportMeasurer
{
public:
	CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept override;
	bool Compare(const CSupport& support1, const CSupport& support2) const noexcept override;
};

// ʹ��MSAC(M-estimator Sample and Consensus)�����е�������Ӧ��������ģ�͵�֧�ֶ�. ���һ��֧�ֶ��и�С��MSAC�÷�, ��ô���֧�ֶȾ͸���
class CMEstimatorSupportMeasurer final : public CSupportMeasurer
{
public:
	CSupport Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept override;
	bool Compare(const CSupport& support1, const CSupport& support2) const noexcept override;
};

