#pragma once
#include <string>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <any>
#include <Windows.h>
#include <math.h>

#include <ceres/ceres.h>

#include <boost/stacktrace.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "Math.h"
#include "Base.h"


struct CBaseOptions
{
	std::string imageDir = "";
	virtual void CheckOptions() const = 0;
};

// SIFT描述子的归一化方式
enum class CSIFTNormalizationType
{
	L1_ROOT, // 先进行L1归一化, 然后求元素的平方根. 这种归一化通常优于标准的L2归一化. 参考文献: "Three things everyone should know to improve object retrieval", Relja Arandjelovic and Andrew Zisserman, CVPR 2012.
	L2  // 传统的L2标准化，即使每个描述子向量的L2范数为1
};
// SIFT特征提取器的选项参数
struct CSIFTExtractionOptions final : public CBaseOptions
{
	// 是否使用GPU
	bool isUseGPU = true;

	// 最大影像尺寸, 如果影像超过这个值, 影像会被缩小
	size_t maxImageSize = 3200;

	// 最多能检测到的特征点数量, 如果检测到的特征点数量超过这个值, 会保留那些尺度更大的特征点
	size_t maxNumFeatures = 8192;

	// 金字塔的总层数
	size_t numOctaves = 4;

	// 每个金字塔层级里的子层个数
	size_t octaveResolution = 3;

	// 用于特征点检测的峰值阈值
	double peakThreshold = 0.02 / octaveResolution;

	// 用于特征点检测的边缘阈值
	double edgeThreshold = 10.0;

	// 是否要对SIFT特征进行仿射形状的估计. 以定向椭圆而不是定向圆盘的形式估算SIFT特征的仿射形状
	bool isEstimateAffineShape = false;

	// 如果isEstimateAffineShape为false的话, 设定每个关键点最多可以有多少个方向
	size_t maxNumOrientations = 2;

	// 固定关键点的方向为0
	bool isUpRight = false;

	// 实现Domain-Size Pooling的相关参数, 计算的是检测的尺度周围多个尺度的平均SIFT描述符
	// 参考文献: "Domain-Size Pooling in Local Descriptors and Network Architectures", J. Dong and S. Soatto, CVPR 2015
	// 参考文献: "Comparative Evaluation of Hand-Crafted and Learned Local Features", Schönberger, Hardmeier, Sattler, Pollefeys, CVPR 2016.
	bool isDomainSizePooling = false;
	double dspMinScale = 1.0 / 6.0;
	double dspMaxScale = 3.0;
	size_t dspNumScales = 10;

	CSIFTNormalizationType normalizationType = CSIFTNormalizationType::L1_ROOT; // SIFT描述子的归一化方式

	inline void CheckOptions() const override
	{
		Check(peakThreshold > 0);
		Check(edgeThreshold > 0);
		Check(dspMinScale > 0);
		Check(dspMaxScale >= dspMinScale);
	}
};

// SIFT特征匹配器的选项参数
struct CSIFTMatchingOptions final : public CBaseOptions
{
	// 是否使用GPU
	bool isUseGPU = true;

	// 最大距离比率: 第一和第二最佳匹配之间的最大距离比. 0.8
	double maxRatio = 0.8;

	// 与最佳匹配之间的最大距离. 值越大, 匹配数越多. 512
	double maxDistance = 128;

	// 是否启用交叉检测. 当启用时, 一个匹配对只有在两个方向(A->B和B->A)都是最佳匹配时, 才算一个有效匹配
	bool doCrossCheck = true;

	// 最大匹配数
	size_t maxNumMatches = 8192;

	// 几何验证时的最大对极误差(单位: 像素). 3
	double maxError = 1;

	// 几何验证的置信度阈值
	double confidence = 0.999;

	// RANSAC迭代的最小和最大次数
	size_t minNumTrials = 100;
	size_t maxNumTrials = 10000;

	// 预先假设的最小内点比率, 用于确定RANSAC的最大迭代次数
	double minInlierRatio = 0.25;

	// 影像对被视作有效的双视几何需要的最小内点数
	size_t minNumInliers = 15;

	// 是否尝试估计多种双视几何模型
	bool isEstimateMultiModels = false;

	// 是否执行引导匹配
	bool isPerformGuidedMatching = false;

	// 是否强制使用单应矩阵来应对平面场景
	bool isPlanarScene = false;

	// 是否估算影像之间的相对位姿并且保存到数据库
	bool isComputeRelativePose = false;

	// 是否使用FLANN算法加速匹配
	bool isUseFLANN = true;

	inline void CheckOptions() const override
	{
		Check(maxRatio > 0);
		Check(maxDistance > 0);
		Check(maxError > 0);
		Check(maxNumTrials > 0);
		Check(maxNumTrials >= minNumTrials);
		Check(minInlierRatio >= 0 && minInlierRatio <= 1);
	}
};

// RANSAC算法的选项参数
struct CRANSACOptions final : public CBaseOptions
{
	// 如果样本的误差的平方不超过maxError, 那么就认为是内点
	double maxError = 4.0;	

	// 预先设定的最小内点比例, 用来计算RANSAC的最大迭代次数. 只会在小于maxNumTrials时适用
	double minInlierRatio = 0.25;

	// 置信度阈值. 如果当前找到一个不含外点的概率达到了这个阈值, 那么终止迭代
	double confidence = 0.999;

	// 基于指定置信值动态计算最大迭代次数的乘数. 基于指定置信度值动态计算的最大迭代次数的num_trials_multiplier的倍数
	double maxIterNumTrialsMultiplier = 3.0;

	// 从随机子集估计模型的最小/最大随机试验次数. RANSAC算法随机选择样本集合进行模型估计的最小/最大次数
	size_t minNumTrials = 100;
	size_t maxNumTrials = 10000;

	inline void CheckOptions() const override
	{
		Check(maxError > 0);
		Check(minInlierRatio >= 0 && minInlierRatio <= 1);
		Check(confidence >= 0 && confidence <= 1);
		Check(minNumTrials <= maxNumTrials);
	}
};

// 双视几何选项参数
struct CTwoViewGeometryOptions final :public CBaseOptions
{
	// 有效双视几何所需要的最小内点数量
	size_t minInliersNum = 15;

	// 通过估计本质矩阵E和基础矩阵F来判断相机是否校正正确
	// 如果E产生的内点数量接近或者超过阈值(=minE_F_InliersRatio*F产生的内点数), 那么就认为相机校正是正确的
	double minE_F_InliersRatio = 0.95;

	// 用于检测退化场景(平面或者纯旋转). 如果单应矩阵H产生的内点数比例接近对极几何的内点数比例(由maxH_InliersRatio控制), 那么就认为是退化的场景
	double maxH_InliersRatio = 0.8;
		
	// 用于检测是否存在水印. 如果"纯平移而且集中在边缘区域的"内点与总内点的比例超过waterMarkMinInlierRatio, 那么就认为存在水印
	double watermarkMinInlierRatio = 0.7;
	
	// 定义图像边界区域的宽度为影像对角线长度的watermarkBorderSize倍
	double watermarkBorderSize = 0.1;

	// 是否启用水印检测. 在图像边界区域(由watermarkBorderSize指定)计算满足纯平移运动条件的内点比率, 如果超过watermarkMinInlierRatio, 就认为存在水印
	bool isDetectWatermark = false;

	// 在多模型估计中是否忽略水印模型
	bool isIgnoreWatermarkInMultiple = true;

	// 是否只估计单应矩阵H(纯平面或者纯旋转场景)
	bool isForceUse_H = false;

	// 是否要计算相对位姿, 会有更多的计算量
	bool isComputeRelativePose = false;

	CRANSACOptions RANSACOptions; // 双视几何估计过程中使用的RANSAC算法的选项参数

	inline void CheckOptions() const override
	{
		Check(minInliersNum >= 0);
		Check(minE_F_InliersRatio >= 0 && minE_F_InliersRatio <= 1);
		Check(maxH_InliersRatio >= 0 && maxH_InliersRatio <= 1);
		Check(watermarkMinInlierRatio >= 0 && watermarkMinInlierRatio <= 1);
		Check(watermarkBorderSize >= 0 && watermarkBorderSize <= 1);
		RANSACOptions.CheckOptions();
	}
};

enum class CTriangulationResidualType
{
	CAngularError,
	CReprojectionError
};
// 三角测量选项参数
struct CTriangulationOptions final :public CBaseOptions
{
	double minTriAngle_Deg = 1.5;
	int maxTransitivity = 1;                  // 寻找对应关系时的最大传递性
	double createMaxAngleError = 2;           // 创建新三角测量时的最大角度误差
	double continueMaxAngleError = 2;         // 继续现有三角测量时的最大角度误差
	double mergeMaxReprojectionError = 4;     // 合并三角测量时的最大像素重投影误差
	double completeMaxReprojectionError = 4;  // 完成现有三角测量时的最大重投影误差
	int completeMaxTransitivity = 5;          // 轨迹完成时的最大传递性
	double retriangulateMaxAngleError = 5.0;  // 重新三角测量时的最大角度误差
	double retriangulateMinRatio = 0.2;       // 重新三角测量时, 影像对之间共同三角测量与对应关系数量的最小比例
	int retriangulateMaxTrials = 1;           // 重新三角测量一个影像对的最大尝试次数
	bool isIgnoreTwoViewTracks = true;        // 是否忽略两视图轨迹

	CTriangulationResidualType residualType = CTriangulationResidualType::CAngularError;
	CRANSACOptions ransacOptions;
	inline void CheckOptions() const override
	{
		Check(minTriAngle_Deg > 0);
		Check(maxTransitivity >= 0);
		Check(createMaxAngleError > 0);
		Check(continueMaxAngleError > 0);
		Check(mergeMaxReprojectionError > 0);
		Check(completeMaxReprojectionError > 0);
		Check(completeMaxTransitivity >= 0);
		Check(retriangulateMaxAngleError > 0);
		Check(retriangulateMinRatio >= 0 && retriangulateMinRatio <= 1);
		Check(retriangulateMaxTrials >= 0);
		ransacOptions.CheckOptions();
	}
};

enum class CNextImageSelectionMethod
{
	CMaxVisiblePointsNum,
	CMaxVisiblePointsRatio,
	CMinUncertainty
};
// 重建选项参数
struct CReconstructionOptions final :public CBaseOptions
{
	size_t minNumMatches = 15;      // 被认为有效的最少匹配数
	size_t maxModelOverlap = 3;     // 子模型之间的最大重叠影像数, 20
	size_t minModelSize = 5;        // 一个子模型所需的最少注册影像数
	size_t numInitialTrials = 20;   // 初始化模型的尝试次数

	// 用于过滤具有退化内参的图像的阈值
	double minFocalLengthRatio = 0.1;
	double maxFocalLengthRatio = 10;
	double maxExtraParam = 1;

	size_t initMinNumInliers = 100;                   // 初始影像对的最小内点数
	double initMaxError = 4;                          // 初始影像对做双视几何估计的最大像素误差
	double initMaxForwardMotion = 0.95;               // 初始影像对的最大前向运动
	double initMinTriAngle = 16;                      // 初始影像对的最小交会角(度)
	size_t initMaxRegTrials = 2;                      // 一张影像最多进行多少次初始化尝试

	double absPoseMaxError = 12;                      // 绝对位姿估计中的最大重投影误差
	size_t absPoseMinNumInliers = 30;                 // 绝对位姿估计中的最小内点数
	double absPoseMinInlierRatio = 0.25;              // 绝对位姿估计中的最小内点比例
	bool isEstimateAbsPoseFocalLength = false;        // 绝对位姿估计中是否估计焦距
	double absPoseMinFocalLengthRatio = 0.1;          // 绝对位姿估计中在给定相机的焦距周围进行离散焦距采样的最小焦距比率
	double absPoseMaxFocalLengthRatio = 10;           // 绝对位姿估计中在给定相机的焦距周围进行离散焦距采样的最大焦距比率
	size_t absPoseNumFocalLengthSamples = 30;         // 绝对位姿估计中用于焦距估计的离散样本数量
	CRANSACOptions absPoseRANSACOoptions;             // 绝对位姿估计中用于P3P RANSAC的选项

	double absPoseRefineGradientTolerance = 1;        // 绝对位姿精化中的收敛准则
	size_t absPoseRefineMaxNumIterations = 100;       // 绝对位姿精化中的求解器的最大迭代次数
	double absPoseRefineLossFunctionScale = 1;        // 绝对位姿精化中的损失函数缩放因子, 决定何时进行残差的鲁棒化
	bool isRefineAbsPoseFocalLength = true;           // 绝对位姿精化中是否优化焦距
	bool isRefineAbsPoseExtraParams = true;           // 绝对位姿精化中是否优化额外相机参数

	size_t numLocalBundleAdjustmentImages = 6;        // 局部平差的影像数量
	double minTriAngleLocalBundleAdjustment = 6;      // 局部平差中选择影像的最小交会角

	double filterMaxReprojectionError = 4;            // 观测值的最大像素重投影误差
	double filterMinTriAngle = 1.5;                   // 稳定3D点的最小交会角
	size_t maxRegTrials = 3;                          // 注册影像的最大尝试次数
	bool isFixingExistingImages = false;              // 如果提供了初始模型, 是否固定其已有影像的姿态

	// 选择下一张最优待注册影像的方法
	CNextImageSelectionMethod nextImageSelectionMethod = CNextImageSelectionMethod::CMinUncertainty;
	inline void CheckOptions() const override
	{
		Check(minNumMatches > 0);
		Check(maxModelOverlap > 0);
		Check(minModelSize > 0);
		Check(numInitialTrials > 0);
		Check(initMinNumInliers > 0);
		Check(initMaxError > 0);
		Check(initMaxForwardMotion >= 0 && initMaxForwardMotion <= 1);
		Check(initMinTriAngle >= 0);
		Check(initMaxRegTrials >= 1);

		Check(absPoseMaxError > 0);
		Check(absPoseMinNumInliers > 0);
		Check(absPoseMinInlierRatio >= 0 && absPoseMinInlierRatio <= 1);
		Check(absPoseNumFocalLengthSamples > 0);
		Check(absPoseMinFocalLengthRatio > 0);
		Check(absPoseMaxFocalLengthRatio > 0);
		Check(absPoseMinFocalLengthRatio < absPoseMaxFocalLengthRatio);
		absPoseRANSACOoptions.CheckOptions();

		Check(absPoseRefineGradientTolerance >= 0);
		Check(absPoseRefineLossFunctionScale >= 0);

		Check(numLocalBundleAdjustmentImages >= 2);
		Check(minTriAngleLocalBundleAdjustment >= 0);
		Check(filterMaxReprojectionError >= 0);
		Check(filterMinTriAngle >= 0);
		Check(maxRegTrials >= 1);
	}
};

enum class CLossFunctionType
{
	CTrivial,  // 非鲁棒的, 对异常值非常敏感
	CSoft_L1,  // 具有一定程度的鲁棒性
	CCauchy    // 鲁棒的, 特别适用于存在许多异常值的情况, 但是下降得相对较慢
};
// 平差选项参数
struct CBundleAdjustmentOptions final :public CBaseOptions
{
	CLossFunctionType lossFunctionType = CLossFunctionType::CTrivial;
	double lossFunctionScale = 1;           // 用于CSoft_L1和CCauchy的缩放因子, 用于决定鲁棒化发生时的残差水平
	bool isRefineFocalLength = true;        // 是否优化焦距
	bool isRefinePrincipalPoint = false;    // 是否优化主点
	bool isRefineExtraParams = true;        // 是否优化相机畸变参数
	bool isRefineExtrinsics = true;         // 是否优化外参
	size_t minNumResidualsForMultiThreading = 50000;  // 启用多线程做平差的最小残差数
	ceres::Solver::Options ceresSolverOptions;        // ceres-solver选项

	CBundleAdjustmentOptions()
	{
		ceresSolverOptions.function_tolerance = 0;
		ceresSolverOptions.gradient_tolerance = 0;
		ceresSolverOptions.parameter_tolerance = 0;
		ceresSolverOptions.minimizer_progress_to_stdout = true;
		ceresSolverOptions.max_num_iterations = 100;
		ceresSolverOptions.max_linear_solver_iterations = 200;
		ceresSolverOptions.max_num_consecutive_invalid_steps = 10;
		ceresSolverOptions.max_consecutive_nonmonotonic_steps = 10;
		ceresSolverOptions.num_threads = -1;
	}

	// 根据指定的选项创建一个新的损失函数. 调用者负责该损失函数的所有权
	inline ceres::LossFunction* CreateLossFunction() const
	{
		ceres::LossFunction* lossFunction = nullptr;
		switch (lossFunctionType)
		{
		case CLossFunctionType::CTrivial:
			lossFunction = new ceres::TrivialLoss();
			break;
		case CLossFunctionType::CSoft_L1:
			lossFunction = new ceres::SoftLOneLoss(lossFunctionScale);
			break;
		case CLossFunctionType::CCauchy:
			lossFunction = new ceres::CauchyLoss(lossFunctionScale);
			break;
		}
		Check(lossFunction);
		return lossFunction;
	}

	inline void CheckOptions() const override
	{
		Check(lossFunctionScale >= 0);
	}
};




struct COptions final : public CBaseOptions
{
	CSIFTExtractionOptions SIFTExtractionOptions;
	CSIFTMatchingOptions SIFTMatchingOptions;
	CRANSACOptions RANSACOptions;
	CTwoViewGeometryOptions twoViewGeometryOptions;
	CTriangulationOptions estimateTriangulationOptions;
	CReconstructionOptions reconstructionOptions;
	CBundleAdjustmentOptions bundleAdjustmentOptions;

	inline void CheckOptions() const override
	{
		SIFTExtractionOptions.CheckOptions();
		SIFTMatchingOptions.CheckOptions();
		RANSACOptions.CheckOptions();
		twoViewGeometryOptions.CheckOptions();
		estimateTriangulationOptions.CheckOptions();
		reconstructionOptions.CheckOptions();
		bundleAdjustmentOptions.CheckOptions();
	}
};















