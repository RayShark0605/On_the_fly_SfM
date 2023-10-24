#include "SupportMeasurer.h"


using namespace std;
CSupport CInlierSupportMeasurer::Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept
{
	CSupport support;
    support.numInliers = 0;
    support.residualSum = 0;
    for (const double residual : residuals) 
    {
        if (residual <= maxResidual)
        {
            support.numInliers++;
            support.residualSum += residual;
        }
    }
    return support;
}
bool CInlierSupportMeasurer::Compare(const CSupport& support1, const CSupport& support2) const noexcept
{
    if (support1.numInliers > support2.numInliers)
    {
        return true;
    }
    return (support1.numInliers == support2.numInliers && support1.residualSum < support2.residualSum);
}

CSupport CMEstimatorSupportMeasurer::Evaluate(const std::vector<double>& residuals, double maxResidual) const noexcept
{
    CSupport support;
    support.numInliers = 0;
    support.score = 0;
    for (const double residual : residuals)
    {
        if (residual <= maxResidual)
        {
            support.numInliers++;
            support.score += residual;
        }
        else
        {
            support.score += maxResidual;
        }
    }
    return support;
}
bool CMEstimatorSupportMeasurer::Compare(const CSupport& support1, const CSupport& support2) const noexcept
{
    return support1.score < support2.score;
}
