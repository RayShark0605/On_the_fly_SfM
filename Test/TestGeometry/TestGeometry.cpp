#include <gtest/gtest.h>
#include <mutex>

#include "../../src/Geometry/Rigid3D.h"
#include "../../src/Geometry/Sim3D.h"
#include "../../src/Geometry/EssentialMatrix.h"
#include "../../src/Geometry/HomographyMatrix.h"
#include "../../src/Geometry/Triangulation.h"

using namespace std;

TEST(DecomposeEssentialMatrix, Nominal)
{
    const CRigid3D camera1ToCamera2(Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d(0.5, 1, 1).normalized());
    const Eigen::Matrix3d camera1ToCamera2_R = camera1ToCamera2.rotation.toRotationMatrix();
    const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2);

    Eigen::Matrix3d R1;
    Eigen::Matrix3d R2;
    Eigen::Vector3d t;
    DecomposeEssentialMatrix(E, R1, R2, t);

    EXPECT_TRUE((R1 - camera1ToCamera2_R).norm() < 1e-10 || (R2 - camera1ToCamera2_R).norm() < 1e-10);
    EXPECT_TRUE((t - camera1ToCamera2.translation).norm() < 1e-10 || (t + camera1ToCamera2.translation).norm() < 1e-10);
}

TEST(PoseToEssentialMatrix, Nominal)
{
    EXPECT_EQ(PoseToEssentialMatrix(CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, 1))), (Eigen::MatrixXd(3, 3) << 0, -1, 0, 1, 0, 0, 0, 0, 0).finished());
    EXPECT_EQ(PoseToEssentialMatrix(CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, 2))), (Eigen::MatrixXd(3, 3) << 0, -1, 0, 1, 0, 0, 0, 0, 0).finished());
}

TEST(PoseFromEssentialMatrix, Nominal)
{
    const CRigid3D worldToCamera1;
    const CRigid3D worldToCamera2(Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 0, 0).normalized());
    const CRigid3D camera1ToCamera2 = worldToCamera2 * worldToCamera1.Inverse();
    const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2);

    vector<Eigen::Vector3d> points3D(4);
    points3D[0] = Eigen::Vector3d(0, 0, 1);
    points3D[1] = Eigen::Vector3d(0, 0.1, 1);
    points3D[2] = Eigen::Vector3d(0.1, 0, 1);
    points3D[3] = Eigen::Vector3d(0.1, 0.1, 1);

    vector<Eigen::Vector2d> points1(4);
    vector<Eigen::Vector2d> points2(4);
    for (size_t i = 0; i < points3D.size(); i++)
    {
        points1[i] = (worldToCamera1 * points3D[i]).hnormalized();
        points2[i] = (worldToCamera2 * points3D[i]).hnormalized();
    }

    points3D.clear();

    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    EssentialMatrixToPose(E, points1, points2, R, t, points3D);

    EXPECT_EQ(points3D.size(), 4);

    EXPECT_TRUE(R.isApprox(camera1ToCamera2.rotation.toRotationMatrix()));
    EXPECT_TRUE(t.isApprox(camera1ToCamera2.translation));
}

TEST(FindOptimalImageObservations, Nominal)
{

    const CRigid3D worldToCamera1;
    const CRigid3D worldToCamera2(Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 0, 0).normalized());
    const Eigen::Matrix3d E = PoseToEssentialMatrix(worldToCamera2 * worldToCamera1.Inverse());

    vector<Eigen::Vector3d> points3D(4);
    points3D[0] = Eigen::Vector3d(0, 0, 1);
    points3D[1] = Eigen::Vector3d(0, 0.1, 1);
    points3D[2] = Eigen::Vector3d(0.1, 0, 1);
    points3D[3] = Eigen::Vector3d(0.1, 0.1, 1);

    // Test if perfect projection is equivalent to optimal image observations.
    for (size_t i = 0; i < points3D.size(); i++)
    {
        const Eigen::Vector2d point1 = (worldToCamera1 * points3D[i]).hnormalized();
        const Eigen::Vector2d point2 = (worldToCamera2 * points3D[i]).hnormalized();
        Eigen::Vector2d optimalPoint1;
        Eigen::Vector2d optimalPoint2;
        FindOptimalImageObservations(E, point1, point2, optimalPoint1, optimalPoint2);
        EXPECT_TRUE(point1.isApprox(optimalPoint1));
        EXPECT_TRUE(point2.isApprox(optimalPoint2));
    }
}

TEST(EpipoleFromEssentialMatrix, Nominal)
{
    const CRigid3D camera1ToCamera2(Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, -1).normalized());
    const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2);

    const Eigen::Vector3d leftEpipole = EssentialMatrixToEpipole(E, true);
    const Eigen::Vector3d rightEpipole = EssentialMatrixToEpipole(E, false);
    EXPECT_TRUE(leftEpipole.isApprox(Eigen::Vector3d(0, 0, 1)));
    EXPECT_TRUE(rightEpipole.isApprox(Eigen::Vector3d(0, 0, 1)));
}

TEST(InvertEssentialMatrix, Nominal)
{
    for (size_t i = 1; i < 10000; i++)
    {
        const CRigid3D camera1ToCamera2(Eigen::Quaterniond(EulerAnglesToRotationMatrix(0, 0.1, 0)), Eigen::Vector3d(0, 0, i).normalized());
        const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2);
        const Eigen::Matrix3d inv_inv_E = InvertEssentialMatrix(InvertEssentialMatrix(E));
        EXPECT_TRUE(E.isApprox(inv_inv_E));
    }
}

TEST(DecomposeHomographyMatrix, Nominal)
{
    Eigen::Matrix3d H;
    H << 2.649157564634028, 4.583875997496426, 70.694447785121326,
        -1.072756858861583, 3.533262150437228, 1513.656999614321649,
        0.001303887589576, 0.003042206876298, 1;
    H *= 3;

    Eigen::Matrix3d K;
    K << 640, 0, 320, 0, 640, 240, 0, 0, 1;

    vector<Eigen::Matrix3d> R;
    vector<Eigen::Vector3d> t;
    vector<Eigen::Vector3d> n;
    DecomposeHomographyMatrix(H, K, K, R, t, n);

    EXPECT_EQ(R.size(), 4);
    EXPECT_EQ(t.size(), 4);
    EXPECT_EQ(n.size(), 4);

    Eigen::Matrix3d R_ref;
    R_ref << 0.43307983549125, 0.545749113549648, -0.717356090899523,
        -0.85630229674426, 0.497582023798831, -0.138414255706431,
        0.281404038139784, 0.67421809131173, 0.682818960388909;
    const Eigen::Vector3d t_ref(1.826751712278038, 1.264718492450820, 0.195080809998819);
    const Eigen::Vector3d n_ref(-0.244875830334816, -0.480857890778889, -0.841909446789566);

    bool isRefSolutionExists = false;
    for (size_t i = 0; i < 4; i++)
    {
        const double kEps = 1e-6;
        if ((R[i] - R_ref).norm() < kEps && (t[i] - t_ref).norm() < kEps && (n[i] - n_ref).norm() < kEps)
        {
            isRefSolutionExists = true;
        }
    }
    EXPECT_TRUE(isRefSolutionExists);
}

TEST(DecomposeHomographyMatrix, Random)
{
    const int numIters = 40000;
    const double epsilon = 1e-6;
    const Eigen::Matrix3d id3 = Eigen::Matrix3d::Identity();
    for (int i = 0; i < numIters; i++)
    {
        const Eigen::Matrix3d H = Eigen::Matrix3d::Random();
        if (abs(H.determinant()) < epsilon)
        {
            continue;
        }
        vector<Eigen::Matrix3d> R;
        vector<Eigen::Vector3d> t;
        vector<Eigen::Vector3d> n;
        DecomposeHomographyMatrix(H, id3, id3, R, t, n);

        EXPECT_EQ(R.size(), 4);
        EXPECT_EQ(t.size(), 4);
        EXPECT_EQ(n.size(), 4);

        // Test that each candidate rotation is a rotation
        for (const Eigen::Matrix3d& candidateR : R)
        {
            const Eigen::Matrix3d orthogError = candidateR.transpose() * candidateR - id3;

            // Check that candidate_R is an orthognal matrix
            EXPECT_LT(orthogError.lpNorm<Eigen::Infinity>(), epsilon);

            // Check determinant is 1
            EXPECT_NEAR(candidateR.determinant(), 1.0, epsilon);
        }
    }
}

TEST(PoseFromHomographyMatrix, Nominal)
{
    const Eigen::Matrix3d K1 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d K2 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d R_ref = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d t_ref(1, 0, 0);
    const Eigen::Vector3d n_ref(-1, 0, 0);
    const double d_ref = 1;
    const Eigen::Matrix3d H = PoseToHomographyMatrix(K1, K2, R_ref, t_ref, n_ref, d_ref);

    vector<Eigen::Vector2d> points1;
    points1.emplace_back(0.1, 0.4);
    points1.emplace_back(0.2, 0.3);
    points1.emplace_back(0.3, 0.2);
    points1.emplace_back(0.4, 0.1);

    vector<Eigen::Vector2d> points2;
    for (const auto& point1 : points1)
    {
        const Eigen::Vector3d point2 = H * point1.homogeneous();
        points2.push_back(point2.hnormalized());
    }

    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Vector3d n;
    vector<Eigen::Vector3d> points3D;
    HomographyMatrixToPose(H, K1, K2, points1, points2, R, t, n, points3D);

    EXPECT_EQ(R, R_ref);
    EXPECT_EQ(t, t_ref);
    EXPECT_EQ(n, n_ref);
    EXPECT_EQ(points3D.size(), points1.size());
}

TEST(HomographyMatrixFromPose, PureRotation)
{
    const Eigen::Matrix3d K1 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d K2 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d t(0, 0, 0);
    const Eigen::Vector3d n(-1, 0, 0);
    const double d = 1;
    const Eigen::Matrix3d H = PoseToHomographyMatrix(K1, K2, R, t, n, d);
    EXPECT_EQ(H, Eigen::Matrix3d::Identity());
}

TEST(HomographyMatrixFromPose, PlanarScene)
{
    const Eigen::Matrix3d K1 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d K2 = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d t(1, 0, 0);
    const Eigen::Vector3d n(-1, 0, 0);
    const double d = 1;
    const Eigen::Matrix3d H = PoseToHomographyMatrix(K1, K2, R, t, n, d);
    Eigen::Matrix3d H_ref;
    H_ref << 2, 0, 0, 0, 1, 0, 0, 0, 1;
    EXPECT_EQ(H, H_ref);
}

TEST(ComputeClosestRotationMatrix, Nominal)
{
    const Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
    EXPECT_LT((ComputeClosestRotationMatrix(A) - A).norm(), 1e-6);
    EXPECT_LT((ComputeClosestRotationMatrix(2 * A) - A).norm(), 1e-6);
}

TEST(DecomposeProjectionMatrix, Nominal)
{
    for (int i = 1; i < 10000; i++)
    {
        Eigen::Matrix3d refK = i * Eigen::Matrix3d::Identity();
        refK(0, 2) = i;
        refK(1, 2) = 2 * i;
        const CRigid3D worldToCamera(Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d::Random());
        const Eigen::Matrix3x4d P = refK * worldToCamera.ToMatrix();
        Eigen::Matrix3d K;
        Eigen::Matrix3d R;
        Eigen::Vector3d T;
        DecomposeProjectionMatrix(P, K, R, T);
        EXPECT_TRUE(refK.isApprox(K, 1e-6));
        EXPECT_TRUE(worldToCamera.rotation.toRotationMatrix().isApprox(R, 1e-6));
        EXPECT_TRUE(worldToCamera.translation.isApprox(T, 1e-6));
    }
}

TEST(CrossProductMatrix, Nominal)
{
    EXPECT_EQ(CrossProductMatrix(Eigen::Vector3d(0, 0, 0)), Eigen::Matrix3d::Zero());
    Eigen::Matrix3d refMatrix;
    refMatrix << 0, -3, 2, 3, 0, -1, -2, 1, 0;
    EXPECT_EQ(CrossProductMatrix(Eigen::Vector3d(1, 2, 3)), refMatrix);
}

TEST(EulerAngles, X)
{
    const double rx = 0.3;
    const double ry = 0;
    const double rz = 0;
    double rxx, ryy, rzz;

    RotationMatrixToEulerAngles(EulerAnglesToRotationMatrix(rx, ry, rz), rxx, ryy, rzz);

    EXPECT_NEAR(rx, rxx, 1e-6);
    EXPECT_NEAR(ry, ryy, 1e-6);
    EXPECT_NEAR(rz, rzz, 1e-6);
}

TEST(EulerAngles, Y)
{
    const double rx = 0;
    const double ry = 0.3;
    const double rz = 0;
    double rxx, ryy, rzz;

    RotationMatrixToEulerAngles(EulerAnglesToRotationMatrix(rx, ry, rz), rxx, ryy, rzz);

    EXPECT_NEAR(rx, rxx, 1e-6);
    EXPECT_NEAR(ry, ryy, 1e-6);
    EXPECT_NEAR(rz, rzz, 1e-6);
}

TEST(EulerAngles, Z)
{
    const double rx = 0;
    const double ry = 0;
    const double rz = 0.3;
    double rxx, ryy, rzz;

    RotationMatrixToEulerAngles(EulerAnglesToRotationMatrix(rx, ry, rz), rxx, ryy, rzz);

    EXPECT_NEAR(rx, rxx, 1e-6);
    EXPECT_NEAR(ry, ryy, 1e-6);
    EXPECT_NEAR(rz, rzz, 1e-6);
}

TEST(EulerAngles, XYZ)
{
    const double rx = 0.1;
    const double ry = 0.2;
    const double rz = 0.3;
    double rxx, ryy, rzz;

    RotationMatrixToEulerAngles(EulerAnglesToRotationMatrix(rx, ry, rz), rxx, ryy, rzz);

    EXPECT_NEAR(rx, rxx, 1e-6);
    EXPECT_NEAR(ry, ryy, 1e-6);
    EXPECT_NEAR(rz, rzz, 1e-6);
}

TEST(AverageQuaternions, Nominal)
{
    vector<Eigen::Quaterniond> quats;
    vector<double> weights;

    quats = { {Eigen::Quaterniond::Identity()} };
    weights = { 1.0 };
    EXPECT_EQ(AverageQuaternions(quats, weights).coeffs(), Eigen::Quaterniond::Identity().coeffs());

    quats = { Eigen::Quaterniond::Identity() };
    weights = { 2.0 };
    EXPECT_EQ(AverageQuaternions(quats, weights).coeffs(), Eigen::Quaterniond::Identity().coeffs());

    quats = { Eigen::Quaterniond::Identity(), Eigen::Quaterniond::Identity() };
    weights = { 1.0, 1.0 };
    EXPECT_EQ(AverageQuaternions(quats, weights).coeffs(), Eigen::Quaterniond::Identity().coeffs());

    quats = { Eigen::Quaterniond::Identity(), Eigen::Quaterniond::Identity() };
    weights = { 1.0, 2.0 };
    EXPECT_EQ(AverageQuaternions(quats, weights).coeffs(), Eigen::Quaterniond::Identity().coeffs());

    quats = { Eigen::Quaterniond::Identity(), Eigen::Quaterniond(2, 0, 0, 0) };
    weights = { 1.0, 2.0 };
    EXPECT_EQ(AverageQuaternions(quats, weights).coeffs(), Eigen::Quaterniond::Identity().coeffs());

    quats = { Eigen::Quaterniond::Identity(), Eigen::Quaterniond(1, 1, 0, 0) };
    weights = { 1.0, 1.0 };
    EXPECT_TRUE(AverageQuaternions(quats, weights).isApprox(Eigen::Quaterniond(0.92388, 0.382683, 0, 0), 1e-6));

    quats = { Eigen::Quaterniond::Identity(), Eigen::Quaterniond(1, 1, 0, 0) };
    weights = { 1.0, 2.0 };
    EXPECT_TRUE(AverageQuaternions(quats, weights).isApprox(Eigen::Quaterniond(0.850651, 0.525731, 0, 0), 1e-6));
}

TEST(InterpolateCameraPoses, Nominal)
{
    const CRigid3D worldToCamera1(Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d::Random());
    const CRigid3D worldToCamera2(Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d::Random());

    const CRigid3D interpolateWorldToCamera1 = InterpolateCameraPoses(worldToCamera1, worldToCamera2, 0);
    EXPECT_TRUE(interpolateWorldToCamera1.translation.isApprox(worldToCamera1.translation));

    const CRigid3D interpolateWorldToCamera2 =InterpolateCameraPoses(worldToCamera1, worldToCamera2, 1);
    EXPECT_TRUE( interpolateWorldToCamera2.translation.isApprox(worldToCamera2.translation));

    const CRigid3D interpolateWorldToCamera3 = InterpolateCameraPoses(worldToCamera1, worldToCamera2, 0.5);
    EXPECT_TRUE(interpolateWorldToCamera3.translation.isApprox( (worldToCamera1.translation + worldToCamera2.translation) / 2));
}

TEST(CheckCheirality, Nominal)
{
    const Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d t(1, 0, 0);

    vector<Eigen::Vector2d> points1;
    vector<Eigen::Vector2d> points2;
    vector<Eigen::Vector3d> points3D;

    points1.emplace_back(0, 0);
    points2.emplace_back(0.1, 0);
    EXPECT_TRUE(CheckCheirality(R, t, points1, points2, points3D));
    EXPECT_EQ(points3D.size(), 1);

    points1.emplace_back(0, 0);
    points2.emplace_back(-0.1, 0);
    EXPECT_TRUE(CheckCheirality(R, t, points1, points2, points3D));
    EXPECT_EQ(points3D.size(), 1);

    points2[1][0] = 0.2;
    EXPECT_TRUE(CheckCheirality(R, t, points1, points2, points3D));
    EXPECT_EQ(points3D.size(), 2);

    points2[0][0] = -0.2;
    points2[1][0] = -0.2;
    EXPECT_FALSE(CheckCheirality(R, t, points1, points2, points3D));
    EXPECT_EQ(points3D.size(), 0);
}

CRigid3D TestRigid3d()
{
    return CRigid3D(Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d::Random());
}
TEST(CRigid3D, Default)
{
    const CRigid3D tform;
    EXPECT_EQ(tform.rotation.coeffs(), Eigen::Quaterniond::Identity().coeffs());
    EXPECT_EQ(tform.translation, Eigen::Vector3d::Zero());
}

TEST(CRigid3D, Inverse)
{
    const CRigid3D aTob = TestRigid3d();
    const CRigid3D bToa = aTob.Inverse();
    for (int i = 0; i < 100; i++)
    {
        const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
        const Eigen::Vector3d x_in_b = aTob * x_in_a;
        EXPECT_LT((bToa * x_in_b - x_in_a).norm(), 1e-6);
    }
}

TEST(CRigid3D, Matrix)
{
    const CRigid3D aTob = TestRigid3d();
    const Eigen::Matrix3x4d aTob_Matrix = aTob.ToMatrix();
    for (int i = 0; i < 100; i++)
    {
        const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
        EXPECT_LT((aTob * x_in_a - aTob_Matrix * x_in_a.homogeneous()).norm(), 1e-6);
    }
}

TEST(CRigid3D, ApplyNoRotation)
{

    const CRigid3D aTob(Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 2, 3));
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(2, 4, 6)).norm(), 1e-6);
}

TEST(CRigid3D, ApplyNoTranslation)
{
    const CRigid3D aTob(Eigen::Quaterniond(Eigen::AngleAxisd( EIGEN_PI / 2, Eigen::Vector3d::UnitX())), Eigen::Vector3d::Zero());
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(1, -3, 2)).norm(), 1e-6);
}

TEST(CRigid3D, ApplyRotationTranslation)
{
    const CRigid3D aTob(Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI / 2, Eigen::Vector3d::UnitX())), Eigen::Vector3d(1, 2, 3));
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(2, -1, 5)).norm(), 1e-6);
}

TEST(CRigid3D, ApplyChain)
{
    const CRigid3D aTob = TestRigid3d();
    const CRigid3D bToc = TestRigid3d();
    const CRigid3D cTod = TestRigid3d();
    const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
    const Eigen::Vector3d x_in_b = aTob * x_in_a;
    const Eigen::Vector3d x_in_c = bToc * x_in_b;
    const Eigen::Vector3d x_in_d = cTod * x_in_c;
    EXPECT_EQ((cTod * (bToc * (aTob * x_in_a))), x_in_d);
}

TEST(CRigid3D, Compose)
{
    const CRigid3D aTob = TestRigid3d();
    const CRigid3D bToc = TestRigid3d();
    const CRigid3D cTod = TestRigid3d();
    const CRigid3D aTod = cTod * bToc * aTob;
    const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
    const Eigen::Vector3d x_in_b = aTob * x_in_a;
    const Eigen::Vector3d x_in_c = bToc * x_in_b;
    const Eigen::Vector3d x_in_d = cTod * x_in_c;
    EXPECT_LT((aTod * x_in_a - x_in_d).norm(), 1e-6);
}

CSim3D TestSim3d()
{
    return CSim3D(GetRandomUniformReal<double>(0.1, 10), Eigen::Quaterniond::UnitRandom(), Eigen::Vector3d::Random());
}

TEST(CSim3D, Default)
{
    const CSim3D tform;
    EXPECT_EQ(tform.scale, 1);
    EXPECT_EQ(tform.rotation.coeffs(), Eigen::Quaterniond::Identity().coeffs());
    EXPECT_EQ(tform.translation, Eigen::Vector3d::Zero());
}

TEST(CSim3D, Inverse)
{
    const CSim3D aTob = TestSim3d();
    const CSim3D bToa = aTob.Inverse();
    for (int i = 0; i < 10000; i++)
    {
        const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
        const Eigen::Vector3d x_in_b = aTob * x_in_a;
        EXPECT_LT((bToa * x_in_b - x_in_a).norm(), 1e-6);
    }
}

TEST(CSim3D, Matrix)
{
    const CSim3D aTob = TestSim3d();
    const Eigen::Matrix3x4d b_from_a_mat = aTob.ToMatrix();
    for (int i = 0; i < 10000; i++)
    {
        const Eigen::Vector3d x_in_a = Eigen::Vector3d::Random();
        EXPECT_LT((aTob * x_in_a - b_from_a_mat * x_in_a.homogeneous()).norm(), 1e-6);
    }
}

TEST(CSim3D, ApplyScaleOnly)
{
    const CSim3D aTob(2, Eigen::Quaterniond::Identity(), Eigen::Vector3d::Zero());
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(2, 4, 6)).norm(), 1e-6);
}

TEST(CSim3D, ApplyTranslationOnly)
{
    const CSim3D aTob(1, Eigen::Quaterniond::Identity(), Eigen::Vector3d(1, 2, 3));
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(2, 4, 6)).norm(), 1e-6);
}

TEST(CSim3D, ApplyRotationOnly)
{
    const CSim3D aTob(1, Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI / 2, Eigen::Vector3d::UnitX())), Eigen::Vector3d::Zero());
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(1, -3, 2)).norm(), 1e-6);
}

TEST(CSim3D, ApplyScaleRotationTranslation)
{
    const CSim3D aTob(2, Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI / 2, Eigen::Vector3d::UnitX())), Eigen::Vector3d(1, 2, 3));
    EXPECT_LT((aTob * Eigen::Vector3d(1, 2, 3) - Eigen::Vector3d(3, -4, 7)).norm(), 1e-6);
}


void TestEstimationWithNumCoords(const size_t numCoords)
{
    const CSim3D srcToTarget_True = TestSim3d();

    vector<Eigen::Vector3d> src;
    vector<Eigen::Vector3d> dst;
    for (size_t i = 0; i < numCoords; i++)
    {
        src.emplace_back(i, i + 2, i * i);
        dst.push_back(srcToTarget_True * src.back());
    }

    CSim3D srcToTarget;
    EXPECT_TRUE(srcToTarget.Estimate(src, dst));
    EXPECT_NEAR(srcToTarget_True.scale, srcToTarget.scale, 1e-6);
    EXPECT_LT(srcToTarget_True.rotation.angularDistance(srcToTarget.rotation), 1e-6);
    EXPECT_LT((srcToTarget_True.translation - srcToTarget.translation).norm(), 1e-6);
}

TEST(CSim3D, EstimateMinimal) { TestEstimationWithNumCoords(3); }

TEST(CSim3D, EstimateOverDetermined) { TestEstimationWithNumCoords(100); }

TEST(CSim3D, EstimateDegenerate)
{
    vector<Eigen::Vector3d> invalidSrcDst(3, Eigen::Vector3d::Zero());
    EXPECT_FALSE(CSim3D().Estimate(invalidSrcDst, invalidSrcDst));
}
string CreateTestDir()
{
    const testing::TestInfo* test_info = testing::UnitTest::GetInstance()->current_test_info();
    ostringstream test_name_stream;
    test_name_stream << test_info->test_suite_name() << "." << test_info->name();
    const string test_name = test_name_stream.str();

    const boost::filesystem::path test_dir = boost::filesystem::path("colmap_test_tmp_test_data") / test_name;

    // Create directory once. Cleanup artifacts from previous test runs.
    static mutex mutex;
    lock_guard<std::mutex> lock(mutex);
    static set<string> existing_test_names;
    if (existing_test_names.count(test_name) == 0)
    {
        if (boost::filesystem::is_directory(test_dir))
        {
            boost::filesystem::remove_all(test_dir);
        }
        boost::filesystem::create_directories(test_dir);
    }
    existing_test_names.insert(test_name);

    return test_dir.string();
}
TEST(CSim3D, ToFromFile)
{
    const string path = CreateTestDir() + "/file.txt";
    const CSim3D written = TestSim3d();
    written.SaveToFile(path);
    CSim3D read;
    read.ReadFromFile(path);
    EXPECT_EQ(written.scale, read.scale);
    EXPECT_EQ(written.rotation.coeffs(), read.rotation.coeffs());
    EXPECT_EQ(written.translation, read.translation);
}

TEST(TriangulatePoint, Nominal)
{
    const vector<Eigen::Vector3d> points3D =
    {
        Eigen::Vector3d(0, 0.1, 0.1),
        Eigen::Vector3d(0, 1, 3),
        Eigen::Vector3d(0, 1, 2),
        Eigen::Vector3d(0.01, 0.2, 3),
        Eigen::Vector3d(-1, 0.1, 1),
        Eigen::Vector3d(0.1, 0.1, 0.2),
    };

    const CRigid3D worldToCamera1;

    for (int z = 0; z < 5; z++)
    {
        const double qz = z / 5.0;
        for (int tx = 0; tx < 100000; tx += 2)
        {
            const CRigid3D worldToCamera2(Eigen::Quaterniond(0.2, 0.3, 0.4, qz), Eigen::Vector3d(tx, 2, 3));
            for (size_t i = 0; i < points3D.size(); i++)
            {
                const Eigen::Vector3d& point3D = points3D[i];
                const Eigen::Vector3d point2D1 = worldToCamera1 * point3D;
                const Eigen::Vector3d point2D2 = worldToCamera2 * point3D;

                const Eigen::Vector3d tri_point3D = TriangulatePoint(worldToCamera1.ToMatrix(), worldToCamera2.ToMatrix(), point2D1.hnormalized(), point2D2.hnormalized());
                if ((point3D - tri_point3D).norm() >= 1e-10)
                {
                    cout << (point3D - tri_point3D).norm() << endl;
                }
                EXPECT_TRUE((point3D - tri_point3D).norm() < 1e-10);
            }
        }
    }
}

TEST(CalculateTriangulationAngle, Nominal)
{
    const Eigen::Vector3d tvec1(0, 0, 0);
    const Eigen::Vector3d tvec2(0, 1, 0);

    EXPECT_NEAR(CalculateTriangulationAngle(tvec1, tvec2, Eigen::Vector3d(0, 0, 100)), 0.009999666687, 1e-8);
    EXPECT_NEAR(CalculateTriangulationAngle(tvec1, tvec2, Eigen::Vector3d(0, 0, 50)), 0.019997333973, 1e-8);
    EXPECT_NEAR(CalculateTriangulationAngles(tvec1, tvec2, { Eigen::Vector3d(0, 0, 50) })[0], 0.019997333973, 1e-8);
}





int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return 0;
}