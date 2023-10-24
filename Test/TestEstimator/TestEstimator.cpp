#include <gtest/gtest.h>

#include "../../src/Estimator/Sampler.h"
#include "../../src/Estimator/SupportMeasurer.h"
#include "../../src/Geometry/Rigid3D.h"
#include "../../src/Estimator/RANSAC.h"
#include "../../src/Geometry/EssentialMatrix.h"
using namespace std;

TEST(CCombinationSampler, LessSamples)
{
    CCombinationSampler sampler(2);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), 10);
    vector<unordered_set<size_t>> sampledSets;
    for (size_t i = 0; i < 10; i++) 
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 2);
        sampledSets.emplace_back(samples.begin(), samples.end());
        EXPECT_EQ(sampledSets.back().size(), 2);
        for (size_t j = 0; j < i; j++) 
        {
            EXPECT_TRUE(sampledSets[j].count(samples[0]) == 0 || sampledSets[j].count(samples[1]) == 0);
        }
    }
    vector<size_t> samples;
    sampler.Sample(samples);
    EXPECT_TRUE(sampledSets[0].count(samples[0]) == 1 && sampledSets[0].count(samples[1]) == 1);
}

TEST(CCombinationSampler, EqualSamples)
{
    CCombinationSampler sampler(5);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), 1);
    for (size_t i = 0; i < 100; i++) 
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 5);
        EXPECT_EQ(unordered_set<size_t>(samples.begin(), samples.end()).size(), 5);
    }
}

TEST(CProgressiveSampler, LessSamples) 
{
    CProgressiveSampler sampler(2);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), numeric_limits<size_t>::max());
    for (size_t i = 0; i < 100; i++) 
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 2);
        EXPECT_EQ(unordered_set<size_t>(samples.begin(), samples.end()).size(), 2);
    }
}

TEST(CProgressiveSampler, EqualSamples)
{
    CProgressiveSampler sampler(5);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), numeric_limits<size_t>::max());
    for (size_t i = 0; i < 100; i++)
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 5);
        EXPECT_EQ(unordered_set<size_t>(samples.begin(), samples.end()).size(),5);
    }
}

TEST(CProgressiveSampler, Progressive)
{
    const size_t kNumSamples = 5;
    CProgressiveSampler sampler(kNumSamples);
    sampler.Initialize(50);
    size_t prev_last_sample = 5;
    for (size_t i = 0; i < 100; i++)
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        for (size_t i = 0; i < samples.size() - 1; i++)
        {
            EXPECT_LT(samples[i], samples.back());
            EXPECT_GE(samples.back(), prev_last_sample);
            prev_last_sample = samples.back();
        }
    }
}

TEST(CRandomSampler, LessSamples)
{
    CRandomSampler sampler(2);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), numeric_limits<size_t>::max());
    for (size_t i = 0; i < 100; i++) 
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 2);
        EXPECT_EQ(unordered_set<size_t>(samples.begin(), samples.end()).size(), 2);
    }
}

TEST(CRandomSampler, EqualSamples)
{
    CRandomSampler sampler(5);
    sampler.Initialize(5);
    EXPECT_EQ(sampler.GetMaxNumSamples(), numeric_limits<size_t>::max());
    for (size_t i = 0; i < 100; i++) 
    {
        vector<size_t> samples;
        sampler.Sample(samples);
        EXPECT_EQ(samples.size(), 5);
        EXPECT_EQ(unordered_set<size_t>(samples.begin(), samples.end()).size(), 5);
    }
}

TEST(CInlierSupportMeasurer, Nominal)
{
    CSupport support1;
    EXPECT_EQ(support1.numInliers, 0);
    EXPECT_EQ(support1.residualSum, numeric_limits<double>::max());
    CInlierSupportMeasurer measurer;
    vector<double> residuals = { -1.0, 0.0, 1.0, 2.0 };
    support1 = measurer.Evaluate(residuals, 1.0);
    EXPECT_EQ(support1.numInliers, 3);
    EXPECT_EQ(support1.residualSum, 0.0);
    CSupport support2;
    support2.numInliers = 2;
    EXPECT_TRUE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.residualSum = support1.residualSum;
    EXPECT_TRUE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.numInliers = support1.numInliers;
    support2.residualSum += 0.01;
    EXPECT_TRUE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.residualSum -= 0.01;
    EXPECT_FALSE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.residualSum -= 0.01;
    EXPECT_FALSE(measurer.Compare(support1, support2));
    EXPECT_TRUE(measurer.Compare(support2, support1));
}

TEST(CMEstimatorSupportMeasurer, Nominal)
{
    CSupport support1;
    EXPECT_EQ(support1.numInliers, 0);
    EXPECT_EQ(support1.score, numeric_limits<double>::max());
    CMEstimatorSupportMeasurer measurer;
    vector<double> residuals = { -1.0, 0.0, 1.0, 2.0 };
    support1 = measurer.Evaluate(residuals, 1.0);
    EXPECT_EQ(support1.numInliers, 3);
    EXPECT_EQ(support1.score, 1.0);
    CSupport support2 = support1;
    EXPECT_FALSE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.numInliers -= 1;
    support2.score += 0.01;
    EXPECT_TRUE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.score -= 0.01;
    EXPECT_FALSE(measurer.Compare(support1, support2));
    EXPECT_FALSE(measurer.Compare(support2, support1));
    support2.score -= 0.01;
    EXPECT_FALSE(measurer.Compare(support1, support2));
    EXPECT_TRUE(measurer.Compare(support2, support1));
}

TEST(CP3PEstimator, P3P)
{
    const vector<Eigen::Vector3d> points3D = 
    {
        Eigen::Vector3d(1, 1, 1),
        Eigen::Vector3d(0, 1, 1),
        Eigen::Vector3d(3, 1.0, 4),
        Eigen::Vector3d(3, 1.1, 4),
        Eigen::Vector3d(3, 1.2, 4),
        Eigen::Vector3d(3, 1.3, 4),
        Eigen::Vector3d(3, 1.4, 4),
        Eigen::Vector3d(2, 1, 7),
    };

    vector<Eigen::Vector3d> points3DFaulty = points3D;
    for (size_t i = 0; i < points3D.size(); i++) 
    {
        points3DFaulty[i](0) = 20;
    }

    for (double qx = 0; qx < 1; qx += 0.2) 
    {
        for (double tx = 0; tx < 1; tx += 0.1) 
        {
            const CRigid3D worldToCamera_True(Eigen::Quaterniond(1, qx, 0, 0).normalized(), Eigen::Vector3d(tx, 0, 0));

            vector<Eigen::Vector2d> points2D(points3D.size());
            for (size_t i = 0; i < points3D.size(); i++) 
            {
                points2D[i] = (worldToCamera_True * points3D[i]).hnormalized();
            }

            CRANSACOptions options;
            options.maxError = 1e-5;
            CP3PEstimator P3PEstimator;
            CRANSAC ransac(options, &P3PEstimator);
            const CRANSACReport report = ransac.Estimate<Eigen::Vector2d, Eigen::Vector3d, Eigen::Matrix3x4d>(points2D, points3D);

            EXPECT_TRUE(report.isSuccess);
            EXPECT_LT((worldToCamera_True.ToMatrix() - report.model).norm(), 1e-5);

            // Test residuals of exact points.
            vector<double> residuals;
            P3PEstimator.Residuals(TypeVec2AnyVec(points2D), TypeVec2AnyVec(points3D), report.model, residuals);
            for (size_t i = 0; i < residuals.size(); i++) 
            {
                EXPECT_TRUE(residuals[i] < 1e-5);
            }

            // Test residuals of faulty points.
            P3PEstimator.Residuals(TypeVec2AnyVec(points2D), TypeVec2AnyVec(points3DFaulty), report.model, residuals);
            for (size_t i = 0; i < residuals.size(); i++) 
            {
                EXPECT_TRUE(residuals[i] > 0.1);
            }
        }
    }
}

TEST(CEPnPEstimator, EPNP) 
{
    const vector<Eigen::Vector3d> points3D = 
    {
        Eigen::Vector3d(1, 1, 1),
        Eigen::Vector3d(0, 1, 1),
        Eigen::Vector3d(3, 1.0, 4),
        Eigen::Vector3d(3, 1.1, 4),
        Eigen::Vector3d(3, 1.2, 4),
        Eigen::Vector3d(3, 1.3, 4),
        Eigen::Vector3d(3, 1.4, 4),
        Eigen::Vector3d(2, 1, 7),
    };

    auto points3D_faulty = points3D;
    for (size_t i = 0; i < points3D.size(); i++) 
    {
        points3D_faulty[i](0) = 20;
    }

    for (double qx = 0; qx < 1; qx += 0.2)
    {
        for (double tx = 0; tx < 1; tx += 0.1)
        {
            const CRigid3D worldToCamera_True(Eigen::Quaterniond(1, qx, 0, 0).normalized(), Eigen::Vector3d(tx, 0, 0));

            // Project points to camera coordinate system.
            vector<Eigen::Vector2d> points2D;
            for (size_t i = 0; i < points3D.size(); i++)
            {
                points2D.push_back((worldToCamera_True * points3D[i]).hnormalized());
            }

            CRANSACOptions options;
            options.maxError = 1e-5;

            CEPnPEstimator EPnPEstimator;
            CRANSAC ransac(options, &EPnPEstimator);
            const CRANSACReport report = ransac.Estimate<Eigen::Vector2d, Eigen::Vector3d, Eigen::Matrix3x4d>(points2D, points3D);

            EXPECT_TRUE(report.isSuccess);
            EXPECT_LT((worldToCamera_True.ToMatrix() - report.model).norm(), 1e-4);


            vector<double> residuals;
            EPnPEstimator.Residuals(TypeVec2AnyVec(points2D), TypeVec2AnyVec(points3D), report.model, residuals);
            for (size_t i = 0; i < residuals.size(); i++)
            {
                EXPECT_TRUE(residuals[i] < 1e-3);
            }


            EPnPEstimator.Residuals(TypeVec2AnyVec(points2D), TypeVec2AnyVec(points3D_faulty), report.model, residuals);
            for (size_t i = 0; i < residuals.size(); i++) 
            {
                EXPECT_TRUE(residuals[i] > 0.1);
            }
        }
    }
}

TEST(CEPnPEstimator, EPNP_BrokenSolveSignCase)
{
    vector<Eigen::Vector2d> points2D;
    points2D.emplace_back(-2.6783007931074532e-01, 5.3457197430746251e-01);
    points2D.emplace_back(-4.2629907287470264e-01, 7.5623350319519789e-01);
    points2D.emplace_back(-1.6767413005963930e-01, -1.3387172544910089e-01);
    points2D.emplace_back(-5.6616329720373559e-02, 2.3621156497739373e-01);
    points2D.emplace_back(-1.7721225948969935e-01, 2.3395366792735982e-02);
    points2D.emplace_back(-5.1836259886632222e-02, -4.4380694271927049e-02);
    points2D.emplace_back(-3.5897765845560037e-01, 1.6252721078589397e-01);
    points2D.emplace_back(2.7057324473684058e-01, -1.4067450104631887e-01);
    points2D.emplace_back(-2.5811166424334520e-01, 8.0167171300227366e-02);
    points2D.emplace_back(2.0239567448222310e-02, -3.2845953375344145e-01);
    points2D.emplace_back(4.2571014715170657e-01, -2.8321173570154773e-01);
    points2D.emplace_back(-5.4597596412987237e-01, 9.1431935871671977e-02);

    vector<Eigen::Vector3d> points3D;
    points3D.emplace_back(4.4276865308679305e+00, -1.3384364366019632e+00, -3.5997423085253892e+00);
    points3D.emplace_back(2.7278555252512309e+00, -3.8152996187231392e-01, -2.6558518399902824e+00);
    points3D.emplace_back(4.8548566083054894e+00, -1.4756197433631739e+00, -6.8274946022490501e-01);
    points3D.emplace_back(3.1523013527998449e+00, -1.3377020437938025e+00, -1.6443269301929087e+00);
    points3D.emplace_back(3.8551679771512073e+00, -1.0557700545885551e+00, -1.1695994508851486e+00);
    points3D.emplace_back(5.9571373150353812e+00, -2.6120646101684555e+00, -1.0841441206050342e+00);
    points3D.emplace_back(6.3287088499358894e+00, -1.1761274755817175e+00, -2.5951879774151583e+00);
    points3D.emplace_back(2.3005305990121250e+00, -1.4019796626800123e+00, -4.4485464455072321e-01);
    points3D.emplace_back(5.9816859934587354e+00, -1.4211814511691452e+00, -2.0285923889293449e+00);
    points3D.emplace_back(5.2543344690665457e+00, -2.3389255564264144e+00, 4.3708173185524052e-01);
    points3D.emplace_back(3.2181599245991688e+00, -2.8906671988445098e+00, 2.6825718150064348e-01);
    points3D.emplace_back(4.4592895306946758e+00, -9.1235241641579902e-03, -1.6555237117970871e+00);

    CEPnPEstimator EPnPEstimator;
    const vector<Eigen::Matrix3x4d> output = AnyVec2TypeVec<Eigen::Matrix3x4d>(EPnPEstimator.Estimate(TypeVec2AnyVec(points2D), TypeVec2AnyVec(points3D)));

    EXPECT_EQ(output.size(), 1);

    double reprojectError = 0.0;
    for (size_t i = 0; i < points3D.size(); i++)
    {
        reprojectError += ((output[0] * points3D[i].homogeneous()).hnormalized() - points2D[i]).norm();
    }
    EXPECT_TRUE(reprojectError < 0.2);
}

TEST(CEssentialMatrixEstimator_5Points, 5Points)
{
    const double points1Raw[] = {
        0.4964, 1.0577, 0.3650,  -0.0919, -0.5412, 0.0159, -0.5239, 0.9467,
        0.3467, 0.5301, 0.2797,  0.0012,  -0.1986, 0.0460, -0.1622, 0.5347,
        0.0796, 0.2379, -0.3946, 0.7969,  0.2,     0.7,    0.6,     0.3 };

    const double points2Raw[] = {
        0.7570, 2.7340, 0.3961,  0.6981, -0.6014, 0.7110, -0.7385, 2.2712,
        0.4177, 1.2132, 0.3052,  0.4835, -0.2171, 0.5057, -0.2059, 1.1583,
        0.0946, 0.7013, -0.6236, 3.0253, 0.5,     0.9,    0.9,     0.2 };

    const size_t numPoints = 12;

    vector<Eigen::Vector2d> points1(numPoints);
    vector<Eigen::Vector2d> points2(numPoints);
    for (size_t i = 0; i < numPoints; i++) {
        points1[i] = Eigen::Vector2d(points1Raw[2 * i], points1Raw[2 * i + 1]);
        points2[i] = Eigen::Vector2d(points2Raw[2 * i], points2Raw[2 * i + 1]);
    }

    // Enforce repeatable tests
    SetPRNGSeed(0);

    CRANSACOptions options;
    options.maxError = 0.02;
    options.confidence = 0.9999;
    options.minInlierRatio = 0.1;

    CEssentialMatrixEstimator_5Points EssentialMatrixEstimator_5Points;
    CRANSAC ransac(options, &EssentialMatrixEstimator_5Points);
    const CRANSACReport report = ransac.Estimate<Eigen::Vector2d, Eigen::Vector2d, Eigen::Matrix3d>(points1, points2);

    vector<double> residuals;
    EssentialMatrixEstimator_5Points.Residuals(TypeVec2AnyVec(points1), TypeVec2AnyVec(points2), report.model, residuals);

    for (size_t i = 0; i < 10; i++) 
    {
        EXPECT_LE(residuals[i], options.maxError * options.maxError);
    }

    EXPECT_FALSE(report.inlierMask[10]);
    EXPECT_FALSE(report.inlierMask[11]);
}

TEST(CEssentialMatrixEstimator_8Points, 8Points)
{
    const double points1_raw[] = { 1.839035, 1.924743, 0.543582, 0.375221, 0.473240, 0.142522, 0.964910, 0.598376, 0.102388, 0.140092, 15.994343, 9.622164, 0.285901, 0.430055, 0.091150, 0.254594};
    const double points2_raw[] = { 1.002114, 1.129644, 1.521742, 1.846002, 1.084332, 0.275134, 0.293328, 0.588992, 0.839509, 0.087290, 1.779735, 1.116857, 0.878616, 0.602447, 0.642616, 1.028681};

    const size_t NumPoints = 8;
    vector<Eigen::Vector2d> points1(NumPoints);
    vector<Eigen::Vector2d> points2(NumPoints);
    for (size_t i = 0; i < NumPoints; i++) 
    {
        points1[i] = Eigen::Vector2d(points1_raw[2 * i], points1_raw[2 * i + 1]);
        points2[i] = Eigen::Vector2d(points2_raw[2 * i], points2_raw[2 * i + 1]);
    }

    CEssentialMatrixEstimator_8Points EssentialMatrixEstimator_8Points;
    const Eigen::Matrix3d E = AnyVec2TypeVec<Eigen::Matrix3d>(EssentialMatrixEstimator_8Points.Estimate(TypeVec2AnyVec(points1), TypeVec2AnyVec(points2)))[0];

    // Reference values.
    EXPECT_TRUE(abs(E(0, 0) - -0.0368602) < 1e-5);
    EXPECT_TRUE(abs(E(0, 1) - 0.265019) < 1e-5);
    EXPECT_TRUE(abs(E(0, 2) - -0.0625948) < 1e-5);
    EXPECT_TRUE(abs(E(1, 0) - -0.299679) < 1e-5);
    EXPECT_TRUE(abs(E(1, 1) - -0.110667) < 1e-5);
    EXPECT_TRUE(abs(E(1, 2) - 0.147114) < 1e-5);
    EXPECT_TRUE(abs(E(2, 0) - 0.169381) < 1e-5);
    EXPECT_TRUE(abs(E(2, 1) - -0.21072) < 1e-5);
    EXPECT_TRUE(abs(E(2, 2) - -0.00401306) < 1e-5);

    // Check that the internal constraint is satisfied (two singular values equal and one zero).
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(E);
    Eigen::Vector3d s = svd.singularValues();
    EXPECT_TRUE(abs(s(0) - s(1)) < 1e-5);
    EXPECT_TRUE(abs(s(2)) < 1e-5);
}

TEST(CFundamentalMatrixEstimator_7Points, 7Points)
{
    const double points1Raw[] = { 0.4964, 1.0577, 0.3650, -0.0919, -0.5412, 0.0159, -0.5239, 0.9467, 0.3467, 0.5301, 0.2797, 0.0012, -0.1986, 0.0460 };

    const double points2Raw[] = { 0.7570, 2.7340, 0.3961, 0.6981, -0.6014, 0.7110, -0.7385, 2.2712, 0.4177, 1.2132, 0.3052, 0.4835, -0.2171, 0.5057 };

    const size_t numPoints = 7;

    vector<Eigen::Vector2d> points1(numPoints);
    vector<Eigen::Vector2d> points2(numPoints);
    for (size_t i = 0; i < numPoints; i++) 
    {
        points1[i] = Eigen::Vector2d(points1Raw[2 * i], points1Raw[2 * i + 1]);
        points2[i] = Eigen::Vector2d(points2Raw[2 * i], points2Raw[2 * i + 1]);
    }

    CFundamentalMatrixEstimator_7Points estimator;
    const Eigen::Matrix3d F = estimator.Estimate(points1, points2)[0];

    // Reference values obtained from Matlab.
    EXPECT_NEAR(F(0, 0), 4.81441976, 1e-6);
    EXPECT_NEAR(F(0, 1), -8.16978909, 1e-6);
    EXPECT_NEAR(F(0, 2), 6.73133404, 1e-6);
    EXPECT_NEAR(F(1, 0), 5.16247992, 1e-6);
    EXPECT_NEAR(F(1, 1), 0.19325606, 1e-6);
    EXPECT_NEAR(F(1, 2), -2.87239381, 1e-6);
    EXPECT_NEAR(F(2, 0), -9.92570126, 1e-6);
    EXPECT_NEAR(F(2, 1), 3.64159554, 1e-6);
    EXPECT_NEAR(F(2, 2), 1., 1e-6);
}

TEST(CFundamentalMatrixEstimator_8Points, 8Points)
{
    const double points1Raw[] = { 1.839035, 1.924743, 0.543582, 0.375221, 0.473240, 0.142522, 0.964910, 0.598376, 0.102388, 0.140092, 15.994343, 9.622164, 0.285901, 0.430055, 0.091150, 0.254594 };

    const double points2Raw[] = { 1.002114, 1.129644, 1.521742, 1.846002, 1.084332, 0.275134, 0.293328, 0.588992, 0.839509, 0.087290, 1.779735, 1.116857, 0.878616, 0.602447, 0.642616, 1.028681};

    const size_t numPoints = 8;
    vector<Eigen::Vector2d> points1(numPoints);
    vector<Eigen::Vector2d> points2(numPoints);
    for (size_t i = 0; i < numPoints; i++) 
    {
        points1[i] = Eigen::Vector2d(points1Raw[2 * i], points1Raw[2 * i + 1]);
        points2[i] = Eigen::Vector2d(points2Raw[2 * i], points2Raw[2 * i + 1]);
    }

    CFundamentalMatrixEstimator_8Points estimator;
    const Eigen::Matrix3d F = estimator.Estimate(points1, points2)[0];

    // Reference values obtained from Matlab.
    EXPECT_TRUE(abs(F(0, 0) - -0.217859) < 1e-5);
    EXPECT_TRUE(abs(F(0, 1) - 0.419282) < 1e-5);
    EXPECT_TRUE(abs(F(0, 2) - -0.0343075) < 1e-5);
    EXPECT_TRUE(abs(F(1, 0) - -0.0717941) < 1e-5);
    EXPECT_TRUE(abs(F(1, 1) - 0.0451643) < 1e-5);
    EXPECT_TRUE(abs(F(1, 2) - 0.0216073) < 1e-5);
    EXPECT_TRUE(abs(F(2, 0) - 0.248062) < 1e-5);
    EXPECT_TRUE(abs(F(2, 1) - -0.429478) < 1e-5);
    EXPECT_TRUE(abs(F(2, 2) - 0.0221019) < 1e-5);
}

TEST(CHomographyMatrixEstimator, HomographyMatrixEstimate)
{
    for (int x = 0; x < 10; x++) 
    {
        Eigen::Matrix3d H0;
        H0 << x, 0.2, 0.3, 30, 0.2, 0.1, 0.3, 20, 1;

        vector<Eigen::Vector2d> src;
        src.emplace_back(x, 0);
        src.emplace_back(1, 0);
        src.emplace_back(2, 1);
        src.emplace_back(10, 30);

        vector<Eigen::Vector2d> dst;

        for (size_t i = 0; i < 4; i++) 
        {
            const Eigen::Vector3d dsth = H0 * src[i].homogeneous();
            dst.push_back(dsth.hnormalized());
        }

        CHomographyMatrixEstimator HomographyMatrixEstimator;
        const vector<Eigen::Matrix3d> models = HomographyMatrixEstimator.Estimate(src, dst);

        vector<double> residuals;
        HomographyMatrixEstimator.Residuals(src, dst, models[0], residuals);

        for (size_t i = 0; i < 4; i++) 
        {
            EXPECT_TRUE(residuals[i] < 1e-6);
        }
    }
}

TEST(CRANSACReport, Report)
{
    CRANSACReport<Eigen::Matrix<double, 3, 4>> report;
    EXPECT_FALSE(report.isSuccess);
    EXPECT_EQ(report.numTrials, 0);
    EXPECT_EQ(report.support.numInliers, 0);
    EXPECT_EQ(report.support.residualSum, numeric_limits<double>::max());
    EXPECT_EQ(report.inlierMask.size(), 0);
}

TEST(CRANSAC, GetNumTrials)
{
    CRANSACOptions options;
    CSimilarityTransformEstimator SimilarityTransformEstimator;
    CRANSAC ransac(options, &SimilarityTransformEstimator);
    EXPECT_EQ(ransac.GetNumTrials(1, 100, 0.99, 1.0), 4605168);
    EXPECT_EQ(ransac.GetNumTrials(10, 100, 0.99, 1.0), 4603);
    EXPECT_EQ(ransac.GetNumTrials(10, 100, 0.999, 1.0), 6905);
    EXPECT_EQ(ransac.GetNumTrials(10, 100, 0.999, 2.0), 13809);
    EXPECT_EQ(ransac.GetNumTrials(100, 100, 0.99, 1.0), 1);
    EXPECT_EQ(ransac.GetNumTrials(100, 100, 0.999, 1.0), 1);
    EXPECT_EQ(ransac.GetNumTrials(100, 100, 0, 1.0), 1);
}

TEST(CAffineTransformEstimator, Nominal)
{
    for (int x = 0; x < 10; x++)
    {
        Eigen::Matrix<double, 2, 3> A;
        A << x / 10.0, 0.2, 0.3, 30, 0.2, 0.1;

        vector<Eigen::Vector2d> src;
        src.emplace_back(x, 0);
        src.emplace_back(1, 0);
        src.emplace_back(2, 1);

        vector<Eigen::Vector2d> dst;
        for (size_t i = 0; i < 3; i++)
        {
            dst.push_back(A * src[i].homogeneous());
        }

        CAffineTransformEstimator estimator;
        const vector<Eigen::Matrix<double, 2, 3>> models = AnyVec2TypeVec<Eigen::Matrix<double, 2, 3>>(estimator.Estimate(TypeVec2AnyVec(src), TypeVec2AnyVec(dst)));

        EXPECT_EQ(models.size(), 1);

        vector<double> residuals;
        estimator.Residuals(TypeVec2AnyVec(src), TypeVec2AnyVec(dst), models[0], residuals);

        EXPECT_EQ(residuals.size(), 3);

        for (size_t i = 0; i < 3; i++)
        {
            EXPECT_LT(residuals[i], 1e-6);
        }
    }
}

TEST(CTranslationTransformEstimator, Estimate)
{
    SetPRNGSeed(0);

    vector<Eigen::Vector2d> src;
    for (size_t i = 0; i < 100; i++) 
    {
        src.emplace_back(GetRandomUniformReal(-1000.0, 1000.0), GetRandomUniformReal(-1000.0, 1000.0));
    }

    Eigen::Vector2d translation(GetRandomUniformReal(-1000.0, 1000.0), GetRandomUniformReal(-1000.0, 1000.0));

    vector<Eigen::Vector2d> dst;
    for (size_t i = 0; i < src.size(); i++)
    {
        dst.push_back(src[i] + translation);
    }

    CTranslationTransformEstimator translationTransformEstimator;
    const Eigen::Vector2d estimatedTranslation = translationTransformEstimator.Estimate(src, dst)[0];

    EXPECT_NEAR(translation(0), estimatedTranslation(0), 1e-6);
    EXPECT_NEAR(translation(1), estimatedTranslation(1), 1e-6);

    vector<double> residuals;
    translationTransformEstimator.Residuals(src, dst, estimatedTranslation, residuals);

    for (size_t i = 0; i < residuals.size(); i++)
    {
        EXPECT_TRUE(residuals[i] < 1e-6);
    }
}

TEST(CenterAndNormalizeImagePoints, Nominal) 
{
    vector<Eigen::Vector2d> points;
    for (size_t i = 0; i < 11; i++)
    {
        points.emplace_back(i, i);
    }

    vector<Eigen::Vector2d> normed_points;
    Eigen::Matrix3d matrix;
    CenterAndNormalizeImagePoints(points, normed_points, matrix);

    EXPECT_EQ(matrix(0, 0), 0.31622776601683794);
    EXPECT_EQ(matrix(1, 1), 0.31622776601683794);
    EXPECT_EQ(matrix(0, 2), -1.5811388300841898);
    EXPECT_EQ(matrix(1, 2), -1.5811388300841898);

    Eigen::Vector2d mean_point(0, 0);
    for (const auto& point : normed_points) 
    {
        mean_point += point;
    }
    EXPECT_LT(abs(mean_point[0]), 1e-6);
    EXPECT_LT(abs(mean_point[1]), 1e-6);
}

TEST(ComputeSquaredSampsonError, Nominal)
{
    vector<Eigen::Vector2d> points1;
    points1.emplace_back(0, 0);
    points1.emplace_back(0, 0);
    points1.emplace_back(0, 0);
    vector<Eigen::Vector2d> points2;
    points2.emplace_back(2, 0);
    points2.emplace_back(2, 1);
    points2.emplace_back(2, 2);

    const Eigen::Matrix3d E = PoseToEssentialMatrix(CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 0, 0)));

    vector<double> residuals;
    ComputeSquaredSampsonError(points1, points2, E, residuals);

    EXPECT_EQ(residuals.size(), 3);
    EXPECT_EQ(residuals[0], 0);
    EXPECT_EQ(residuals[1], 0.5);
    EXPECT_EQ(residuals[2], 2);
}

TEST(ComputeSquaredReprojectionError, Nominal)
{
    vector<Eigen::Vector2d> points2D;
    points2D.emplace_back(0, 0);
    points2D.emplace_back(0, 0);
    points2D.emplace_back(0, 0);
    vector<Eigen::Vector3d> points3D;
    points3D.emplace_back(2, 0, 1);
    points3D.emplace_back(2, 1, 1);
    points3D.emplace_back(2, 2, 1);

    const CRigid3D cam_from_world(Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 0, 0));

    vector<double> residuals;
    ComputeSquaredReprojectionError(points2D, points3D, cam_from_world.ToMatrix(), residuals);

    EXPECT_EQ(residuals.size(), 3);
    EXPECT_EQ(residuals[0], 9);
    EXPECT_EQ(residuals[1], 10);
    EXPECT_EQ(residuals[2], 13);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return 0;
}