#pragma once
#include <vector>
#include <array>
#include <any>


#include "../Scene/Camera.h"
#include "../Scene/Projection.h"
#include "../Base/Base.h"
#include "../Base/Options.h"
#include "Common.h"

// �����. ������
class CEstimator
{
public:
	const size_t minNumSamples; // ����ģ���������������������
	virtual std::vector<std::any> Estimate(const std::vector<std::any>& X, const std::vector<std::any>& Y) = 0;
	virtual void Residuals(const std::vector<std::any>& X, const std::vector<std::any>& Y, const std::any& model, std::vector<double>& residuals) = 0;

protected:
	explicit CEstimator(size_t minNumSamples);
	size_t GetNumTrials(size_t minNumSamples, size_t numInliers, size_t numSamples, double confidence, double numTrialsMultiplier) const;
};


struct CP3PEstimateRANSACReport
{
	bool isSuccess = false;          // RANSAC�����Ƿ�ɹ�
	size_t numTrials = 0;            // RANSAC�����Ĵ���
	CSupport support;                // ����ģ�͵�֧�ֶ�
	std::vector<char> inlierMask;    // ����, ���ڵ��������Ϊtrue
	Eigen::Matrix3x4d model;         // RANSAC���Ƴ���ģ��
};
// P3P���������. �ο�����: X.S. Gao, X.-R. Hou, J. Tang, H.-F. Chang. Complete Solution Classification for the Perspective-Three-Point Problem.
class CP3PEstimator final :public CEstimator
{
public:
	CP3PEstimator();
	
	// ��3��2D-3D���Ӧ��ϵ�й���P3P����Ľ�
	//
	// @param points2D   ��һ����2Dͼ���. ʵ������Ӧ����vector<Eigen::Vector2d>, vector�Ĵ�СӦ��Ϊ3
	// @param points3D   3D��. ʵ������Ӧ����vector<Eigen::Vector3d>, vector�Ĵ�СӦ��Ϊ3
	//
	// @return           ���ܵ������̬. ʵ������Ӧ����vector<Eigen::Matrix3x4d>
	std::vector<std::any> Estimate(const std::vector<std::any>& points2D, const std::vector<std::any>& points3D) override;
	std::vector<Eigen::Matrix3x4d> Estimate(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D) const;
	CP3PEstimateRANSACReport EstimateRANSAC(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const CRANSACOptions& options, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr) const;
	CP3PEstimateRANSACReport EstimateLoRANSAC(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const CRANSACOptions& options, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr) const;

	// ����һ��2D-3D���Ӧ��ϵ��һ��ͶӰ����, ����ƽ����ͶӰ���
	// 
	// @param points2D       ��һ����2Dͼ���. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points3D       3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	// @param projectMatrix  3��4��ͶӰ����. ʵ������Ӧ����Eigen::Matrix3x4d
	// @param residuals      ���: ƽ����ͶӰ�������
	void Residuals(const std::vector<std::any>& points2D, const std::vector<std::any>& points3D, const std::any& projectMatrix, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& projectMatrix, std::vector<double>& residuals) const;
};

// ʹ��EPnP�������PnP����, ������Ҫ4��2D-3D��Ӧ��ϵ. �ο�����: Lepetit, Vincent, Francesc Moreno-Noguer, and Pascal Fua. "Epnp: An accurate o (n) solution to the pnp problem."
// ��ʵ�ֻ�������ԭʼ�Ŀ�Դ����, ������ֲ��Eigen��������ԭʼ����Ķ���Ľ�
class CEPnPEstimator final : public CEstimator
{
public:
	CEPnPEstimator();

	// �Ӷ��2D-3D���Ӧ��ϵ�й���P3P����Ľ�
	//
	// @param points2D   ��һ����2Dͼ���. ʵ������Ӧ����vector<Eigen::Vector2d>, vector�Ĵ�СӦ��Ϊ3
	// @param points3D   3D��. ʵ������Ӧ����vector<Eigen::Vector3d>, vector�Ĵ�СӦ��Ϊ3
	//
	// @return           ���ܵ������̬. ʵ������Ӧ����vector<Eigen::Matrix3x4d>, vector�Ĵ�СӦ��Ϊ0(���Ʋ��ɹ�)��1
	std::vector<std::any> Estimate(const std::vector<std::any>& points2D, const std::vector<std::any>& points3D) override;
	std::vector<Eigen::Matrix3x4d> Estimate(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D);

	// ����һ��2D-3D���Ӧ��ϵ��һ��ͶӰ����, ����ƽ����ͶӰ���
	// 
	// @param points2D       ��һ����2Dͼ���. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points3D       3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	// @param projectMatrix  3��4��ͶӰ����. ʵ������Ӧ����Eigen::Matrix3x4d
	// @param residuals      ���: ƽ����ͶӰ�������
	void Residuals(const std::vector<std::any>& points2D, const std::vector<std::any>& points3D, const std::any& projectMatrix, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& projectMatrix, std::vector<double>& residuals);

private:
	const std::vector<Eigen::Vector2d>* points2D = nullptr;
	const std::vector<Eigen::Vector3d>* points3D = nullptr;
	std::vector<Eigen::Vector3d> pcs;
	std::vector<Eigen::Vector4d> alphas;
	std::array<Eigen::Vector3d, 4> cws;
	std::array<Eigen::Vector3d, 4> ccs;

	bool ComputePose(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, Eigen::Matrix3x4d& projectMatrix);
	void ChooseControlPoints();
	bool ComputeBarycentricCoordinates();
	Eigen::Matrix<double, Eigen::Dynamic, 12> ComputeM();
	Eigen::Matrix<double, 6, 10> ComputeL6x10(const Eigen::Matrix<double, 12, 12>& Ut);
	Eigen::Matrix<double, 6, 1> ComputeRho();

	// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betasApprox_1  = [B11 B12     B13         B14            ]
	void FindBetasApprox1(const Eigen::Matrix<double, 6, 10>& L_6x10, const Eigen::Matrix<double, 6, 1>& rho, Eigen::Vector4d& betas);

	// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betasApprox_2  = [B11 B12 B22                            ]
	void FindBetasApprox2(const Eigen::Matrix<double, 6, 10>& L_6x10, const Eigen::Matrix<double, 6, 1>& rho, Eigen::Vector4d& betas);

	// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betasApprox_3  = [B11 B12 B22 B13 B23                    ]
	void FindBetasApprox3(const Eigen::Matrix<double, 6, 10>& L_6x10, const Eigen::Matrix<double, 6, 1>& rho, Eigen::Vector4d& betas);
	void RunGaussNewton(const Eigen::Matrix<double, 6, 10>& L_6x10, const Eigen::Matrix<double, 6, 1>& rho, Eigen::Vector4d& betas);
	double ComputeRT(const Eigen::Matrix<double, 12, 12>& Ut, const Eigen::Vector4d& betas, Eigen::Matrix3d& R, Eigen::Vector3d& t);
	void ComputeCcs(const Eigen::Vector4d& betas, const Eigen::Matrix<double, 12, 12>& Ut);
	void ComputePcs();
	void SolveForSign();
	void EstimateRT(Eigen::Matrix3d& R, Eigen::Vector3d& t);
	double ComputeTotalReprojectionError(const Eigen::Matrix3d& R, const Eigen::Vector3d& t);
};

struct CEssentialMatrixEstimate_5PointsRANSACReport
{
	bool isSuccess = false;          // RANSAC�����Ƿ�ɹ�
	size_t numTrials = 0;            // RANSAC�����Ĵ���
	CSupport support;                // ����ģ�͵�֧�ֶ�
	std::vector<char> inlierMask;    // ����, ���ڵ��������Ϊtrue
	Eigen::Matrix3d model;           // RANSAC���Ƴ���ģ��
};
// ʹ����㷨�ӹ�һ��2D��Թ��Ʊ��ʾ���E. �ο�����: D. Nister, An efficient solution to the five-point relative pose problem, IEEE-T-PAMI, 26(6), 2004.
class CEssentialMatrixEstimator_5Points final :public CEstimator
{
public:
	CEssentialMatrixEstimator_5Points();

	// ������5���Ӧ���й���10�����ܵı��ʾ���
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          ���п��ܵ�3��3���ʾ���. ʵ������Ӧ����vector<Eigen::Matrix3d>, vector�Ĵ�СӦ��Ϊ10
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix3d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;
	CEssentialMatrixEstimate_5PointsRANSACReport EstimateLoRANSAC(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const CRANSACOptions& options, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr) const;

	// ����һ���Ӧ��ͱ��ʾ���, ����ƽ����Sampson���
	// 
	// @param points1       ��һ���Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       �ڶ����Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param E             3��3�ı��ʾ���. ʵ������Ӧ����Eigen::Matrix3d
	// @param residuals     ���: ƽ����Sampson���
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& E, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, std::vector<double>& residuals) const;

private:
	void CalculateA(const Eigen::Matrix<double, 9, 4>& E, Eigen::Matrix<double, 10, 20>& A) const;
	void CalculateCoeffs(const Eigen::Matrix<double, 13, 3>& B, Eigen::Matrix<double, 11, 1>& coeffs) const;

};

// ʹ�ð˵㷨�ӹ�һ��2D��Թ��Ʊ��ʾ���E. �ο�����: Hartley and Zisserman, Multiple View Geometry, algorithm 11.1, page 282.
class CEssentialMatrixEstimator_8Points final :public CEstimator
{
public:
	CEssentialMatrixEstimator_8Points();

	// ������8���Ӧ���й��㱾�ʾ���
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          3��3���ʾ���. ʵ������Ӧ����vector<Eigen::Matrix3d>, vector�Ĵ�СӦ��Ϊ1
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix3d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;

	// ����һ���Ӧ��ͱ��ʾ���, ����ƽ����Sampson���
	// 
	// @param points1       ��һ���Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       �ڶ����Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param E             3��3�ı��ʾ���. ʵ������Ӧ����Eigen::Matrix3d
	// @param residuals     ���: ƽ����Sampson���
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& E, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, std::vector<double>& residuals) const;
};

struct CFundamentalMatrixEstimate_7PointsRANSACReport
{
	bool isSuccess = false;          // RANSAC�����Ƿ�ɹ�
	size_t numTrials = 0;            // RANSAC�����Ĵ���
	CSupport support;                // ����ģ�͵�֧�ֶ�
	std::vector<char> inlierMask;    // ����, ���ڵ��������Ϊtrue
	Eigen::Matrix3d model;           // RANSAC���Ƴ���ģ��
};
// ʹ���ߵ㷨���ƻ�������F. �ο�����: Zhengyou Zhang and T. Kanade, Determining the Epipolar Geometry and its Uncertainty: A Review, International Journal of Computer Vision, 1998.
class CFundamentalMatrixEstimator_7Points final :public CEstimator
{
public:
	CFundamentalMatrixEstimator_7Points();

	// ��7�Ե��й���1����3�����ܵĻ�������F
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>, vector�Ĵ�СӦ��Ϊ7
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>, vector�Ĵ�СӦ��Ϊ7
	//
	// @return          3��3��������. ʵ������Ӧ����vector<Eigen::Matrix3d>, vector�Ĵ�С������3
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix3d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;
	CFundamentalMatrixEstimate_7PointsRANSACReport EstimateLoRANSAC(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const CRANSACOptions& options, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr) const;

	// ����һ���Ӧ��ͻ�������, ����ƽ����Sampson���
	// 
	// @param points1       ��һ���Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       �ڶ����Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param F             3��3�Ļ�������. ʵ������Ӧ����Eigen::Matrix3d
	// @param residuals     ���: ƽ����Sampson���
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& F, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& F, std::vector<double>& residuals) const;
};

// ʹ�ð˵㷨���ƻ�������F. �ο�����: Hartley and Zisserman, Multiple View Geometry, algorithm 11.1, page 282.
class CFundamentalMatrixEstimator_8Points final :public CEstimator
{
public:
	CFundamentalMatrixEstimator_8Points();

	// ������8�Ե��й��ƻ�������F
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          3��3��������. ʵ������Ӧ����vector<Eigen::Matrix3d>, vector�Ĵ�СΪ1
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix3d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;

	// ����һ���Ӧ��ͻ�������, ����ƽ����Sampson���
	// 
	// @param points1       ��һ���Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       �ڶ����Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param F             3��3�Ļ�������. ʵ������Ӧ����Eigen::Matrix3d
	// @param residuals     ���: ƽ����Sampson���
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& F, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& F, std::vector<double>& residuals) const;
};

struct CHomographyMatrixEstimateRANSACReport
{
	bool isSuccess = false;          // RANSAC�����Ƿ�ɹ�
	size_t numTrials = 0;            // RANSAC�����Ĵ���
	CSupport support;                // ����ģ�͵�֧�ֶ�
	std::vector<char> inlierMask;    // ����, ���ڵ��������Ϊtrue
	Eigen::Matrix3d model;           // RANSAC���Ƴ���ģ��
};
// ʹ��DLT(ֱ�����Ա任)�㷨���Ƶ�Ӧ����H, ������Ҫ4�Ե���������С���˹���
class CHomographyMatrixEstimator final :public CEstimator
{
public:
	CHomographyMatrixEstimator();

	// ������4�Ե��й��Ƶ�Ӧ����H
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          3��3��α任����. ʵ������Ӧ����vector<Eigen::Matrix3d>, vector�Ĵ�СΪ1
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix3d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;
	CHomographyMatrixEstimateRANSACReport EstimateLoRANSAC(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const CRANSACOptions& options, CSupportMeasurer* supportMeasurer = nullptr, CSampler* sampler = nullptr) const;

	// ����һ���Ӧ��͵�Ӧ����, ����任����ƽ��
	// 
	// @param points1       ��һ���Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       �ڶ����Ӧ��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param H             3��3��ͶӰ����. ʵ������Ӧ����Eigen::Matrix3d
	// @param residuals     ���: �任����ƽ��
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& H, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& H, std::vector<double>& residuals) const;
};

// ���Ʒ���任, ������Ҫ3��2D��
class CAffineTransformEstimator final : public CEstimator
{
public:
	CAffineTransformEstimator();

	// ������3���Ӧ���й��Ʒ���任
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          ���Ƶķ���任����. ʵ������Ӧ����vector<Eigen::Matrix<double, 2, 3>>, vector�Ĵ�СӦ��Ϊ1
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Matrix<double, 2, 3>> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;

	// ����ÿ����ı任����ƽ��
	// 
	// @param points1       2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param A             ����任����. ʵ������Ӧ����Eigen::Matrix<double, 2, 3>
	// @param residuals     ���: ÿ����ı任����ƽ��
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& A, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix<double, 2, 3>& A, std::vector<double>& residuals) const;
};

// ����Դ��Ŀ������ϵͳ�Ķ�Ӧ��, ������ά�����Ա任
// �ο�����: S. Umeyama. Least-Squares Estimation of Transformation Parameters Between Two Point Patterns. IEEE Transactions on Pattern Analysis and Machine Intelligence, Volume 13 Issue 4, Page 376-380, 1991.
class CSimilarityTransformEstimator final : public CEstimator
{
public:
	// isEstimateScale��ʾ"�Ƿ���Ҫ���Ƴ߶�"
	explicit CSimilarityTransformEstimator(bool isEstimateScale = true);

	// ������ά�����Ա任����
	//
	// @param src   3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	// @param dst   3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	//
	// @return      ���Ƶ������Ա任����. ʵ������Ӧ����vector<Eigen::Matrix<double, 3, 4>>, vector�Ĵ�СӦ��Ϊ1
	std::vector<std::any> Estimate(const std::vector<std::any>& src, const std::vector<std::any>& dst) override;
	std::vector<Eigen::Matrix<double, 3, 4>> Estimate(const std::vector<Eigen::Vector3d>& src, const std::vector<Eigen::Vector3d>& dst) const;

	// ����ÿ����ı任�в�(��Դ����ϵ�任��Ŀ������ϵʱ��ƽ���任���)
	// 
	// @param src           3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	// @param dst           3D��. ʵ������Ӧ����vector<Eigen::Vector3d>
	// @param matrix        ����任����. ʵ������Ӧ����Eigen::Matrix<double, 3, 4>
	// @param residuals     ���: ÿ����ı任�в�
	void Residuals(const std::vector<std::any>& src, const std::vector<std::any>& dst, const std::any& matrix, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector3d>& src, const std::vector<Eigen::Vector3d>& dst, const Eigen::Matrix<double, 3, 4>& matrix, std::vector<double>& residuals) const;

private:
	bool isEstimateScale;
};

// ���ƶ�άƽ�Ʊ任, ������Ҫ1��2D��
class CTranslationTransformEstimator final : public CEstimator
{
public:
	CTranslationTransformEstimator();

	// ���ƶ�άƽ�Ʊ任
	//
	// @param points1   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2   2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	//
	// @return          ƽ������. ʵ������Ӧ����vector<Eigen::Vector2d>, vector�Ĵ�СӦ��Ϊ1
	std::vector<std::any> Estimate(const std::vector<std::any>& points1, const std::vector<std::any>& points2) override;
	std::vector<Eigen::Vector2d> Estimate(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2) const;

	// ����ƽ������ƽ��
	// 
	// @param points1       2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param points2       2D��. ʵ������Ӧ����vector<Eigen::Vector2d>
	// @param translation   ƽ������. ʵ������Ӧ����Eigen::Vector2d
	// @param residuals     ���: ÿһ�Ե�Ĳв�
	void Residuals(const std::vector<std::any>& points1, const std::vector<std::any>& points2, const std::any& translation, std::vector<double>& residuals) override;
	void Residuals(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Vector2d& translation, std::vector<double>& residuals) const;
};


struct CTriangulationPoint
{
	Eigen::Vector2d point;            // Ӱ�����ֵ, ������Ϊ��λ. ֻ��Ҫ�ڼ�����ͶӰ���(CReprojectionError)ʱ����
	Eigen::Vector2d pointNormalized;  // ��һ��Ӱ�����ֵ. ʼ�ն�Ҫ����

	CTriangulationPoint() = default;
	CTriangulationPoint(const Eigen::Vector2d& point, const Eigen::Vector2d& pointNormalized)
	{
		this->point = point;
		this->pointNormalized = pointNormalized;
	}
};
struct CTriangulationPose
{
	Eigen::Matrix3x4d projectionMatrix; // �۲�Ӱ���ͶӰ����
	Eigen::Vector3d projectionCenter;   // �۲�Ӱ���ͶӰ����
	const CCamera* camera = nullptr;    // �۲�Ӱ������

	CTriangulationPose() = default;
	CTriangulationPose(const Eigen::Matrix3x4d& projectionMatrix, const Eigen::Vector3d& pose, const CCamera* camera)
	{
		this->projectionMatrix = projectionMatrix;
		projectionCenter = pose;
		this->camera = camera;
	}
};
// ���ǲ������������, �Ӷ���۲�ֵ����3D��, �۲�ֵ��Ӱ�����ֵ�Լ���Ӧ�����λ�˼��ڲ����
// ��������2��Լ������: 1. ����Ǳ����㹻��. 2. ���й۲ⶼ��������������Լ��(cheirality constraint)
class CTriangulationEstimator final : public CEstimator
{
public:
	CTriangulationEstimator(double minTriAngle, CTriangulationResidualType residualType);

	// ���ƶ�άƽ�Ʊ任
	//
	// @param points   Ӱ��۲�ֵ. ʵ������Ӧ����vector<CTriangulationPoint>
	// @param poses    ���λ��. ʵ������Ӧ����vector<CTriangulationPose>
	//
	// @return         3D��Ĺ���ֵ. ʵ������Ӧ����vector<Eigen::Vector3d>, ����ɹ���vector�Ĵ�СӦ��Ϊ1, ʧ����vector�Ĵ�СΪ0
	std::vector<std::any> Estimate(const std::vector<std::any>& points, const std::vector<std::any>& poses) override;
	std::vector<Eigen::Vector3d> Estimate(const std::vector<CTriangulationPoint>& points, const std::vector<CTriangulationPose>& poses) const;

	// ��ƽ����ͶӰ����Ƕ�������ʽ����в�
	// 
	// @param points      Ӱ��۲�ֵ. ʵ������Ӧ����vector<CTriangulationPoint>
	// @param poses       ���λ��. ʵ������Ӧ����vector<CTriangulationPose>
	// @param XYZ         3D��. ʵ������Ӧ����Eigen::Vector3d
	// @param residuals   ���: ÿ�ι۲�Ĳв�
	void Residuals(const std::vector<std::any>& points, const std::vector<std::any>& poses, const std::any& XYZ, std::vector<double>& residuals) override;
	void Residuals(const std::vector<CTriangulationPoint>& points, const std::vector<CTriangulationPose>& poses, const Eigen::Vector3d& XYZ, std::vector<double>& residuals) const;

private:
	CTriangulationResidualType residualType = CTriangulationResidualType::CReprojectionError;
	double minTriAngle = 0.0;
};





