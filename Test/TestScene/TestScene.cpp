#include <limits>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Householder>
#include <Eigen/MetisSupport>
#include <Eigen/Sparse>
#include <Eigen/StdDeque>
#include <Eigen/StdVector>
#include <Eigen/StdList>
#include <Eigen/Geometry>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>


#include "../../src/Scene/Camera.h"
#include "../../src/Scene/CameraModel.h"
#include "../../src/Scene/Image.h"
#include "../../src/Scene/Point3D.h"
#include "../../src/Scene/Projection.h"
#include "../../src/Scene/TwoViewGeometry.h"
#include "../../src/Scene/Database.h"
using namespace std;

TEST(CCamera, Empty)
{
    CCamera camera;
    EXPECT_EQ(camera.GetCameraID(), numeric_limits<size_t>::max());
    EXPECT_EQ(camera.GetCameraModelID(), numeric_limits<size_t>::max());
    EXPECT_EQ(camera.GetCameraModelName(), "Unknown");
    EXPECT_EQ(camera.GetWidth(), 0);
    EXPECT_EQ(camera.GetHeight(), 0);
    EXPECT_FALSE(camera.IsFocalLengthPrior());
    EXPECT_EQ(camera.GetParamsString(), "");
    EXPECT_EQ(camera.GetParams().size(), 0);
    EXPECT_EQ(camera.GetParamsData(), camera.GetParams().data());
}

TEST(CCamera, CameraId)
{
    CCamera camera;
    EXPECT_EQ(camera.GetCameraID(), numeric_limits<size_t>::max());
    camera.SetCameraID(1);
    EXPECT_EQ(camera.GetCameraID(), 1);
}

TEST(CCamera, ModelId)
{
    CCamera camera;
    EXPECT_EQ(camera.GetCameraModelID(), numeric_limits<size_t>::max());
    EXPECT_EQ(camera.GetCameraModelName(), "Unknown");
    camera.SetCameraModelName("Simple Radial");
    EXPECT_EQ(camera.GetCameraModelID(), static_cast<int>((size_t)CCameraModelType::CSimpleRadialCameraModel));
    EXPECT_EQ(camera.GetCameraModelName(), "Simple Radial");
    EXPECT_EQ(camera.GetParamsNum(), 4);
}

TEST(CCamera, WidthHeight)
{
    CCamera camera;
    EXPECT_EQ(camera.GetWidth(), 0);
    EXPECT_EQ(camera.GetHeight(), 0);
    camera.SetWidth(1);
    EXPECT_EQ(camera.GetWidth(), 1);
    EXPECT_EQ(camera.GetHeight(), 0);
    camera.SetHeight(1);
    EXPECT_EQ(camera.GetWidth(), 1);
    EXPECT_EQ(camera.GetHeight(), 1);
}

TEST(CCamera, FocalLength)
{
    CCamera camera;
    camera = CCamera((size_t)CCameraModelType::CSimpleRadialCameraModel, 1.0, 1, 1);
    EXPECT_EQ(camera.GetFocalLength(), 1.0);
    camera.SetFocalLength(2.0);
    EXPECT_EQ(camera.GetFocalLength(), 2.0);
}

TEST(CCamera, PriorFocalLength)
{
    CCamera camera;
    EXPECT_FALSE(camera.IsFocalLengthPrior());
    camera.SetFocalLengthPrior(true);
    EXPECT_TRUE(camera.IsFocalLengthPrior());
    camera.SetFocalLengthPrior(false);
    EXPECT_FALSE(camera.IsFocalLengthPrior());
}

TEST(CCamera, PrincipalPoint)
{
    CCamera camera;
    camera= CCamera((size_t)CCameraModelType::CSimpleRadialCameraModel, 1.0, 1, 1);
    EXPECT_EQ(camera.GetPrincipalPointX(), 0.5);
    EXPECT_EQ(camera.GetPrincipalPointY(), 0.5);
    camera.SetPrincipalPointX(2.0);
    EXPECT_EQ(camera.GetPrincipalPointX(), 2.0);
    EXPECT_EQ(camera.GetPrincipalPointY(), 0.5);
    camera.SetPrincipalPointY(2.0);
    EXPECT_EQ(camera.GetPrincipalPointX(), 2.0);
    EXPECT_EQ(camera.GetPrincipalPointY(), 2.0);
}

TEST(CCamera, ParamsInfo)
{
    CCamera camera;
    camera.SetCameraModelID((size_t)CCameraModelType::CSimpleRadialCameraModel);
    EXPECT_EQ(camera.GetParamsType(), "f, cx, cy, k");
}

TEST(CCamera, IsUndistorted)
{
    CCamera camera;
    camera= CCamera((size_t)CCameraModelType::CSimpleRadialCameraModel, 1.0, 1, 1);
    EXPECT_TRUE(camera.IsUndistorted());
    camera.SetParams({ 1.0, 0.5, 0.5, 0.005 });
    EXPECT_FALSE(camera.IsUndistorted());
}

TEST(CCamera, IsBogusParams)
{
    CCamera camera;
    camera= CCamera((size_t)CCameraModelType::CSimpleRadialCameraModel, 1.0, 1, 1);
    EXPECT_FALSE(camera.IsBogusParams(0.1, 1.1, 1.0));
    camera.GetParams(3) = 1.01;
    EXPECT_TRUE(camera.IsBogusParams(0.1, 1.1, 1.0));
    camera.GetParams(3) = -0.5;
    EXPECT_FALSE(camera.IsBogusParams(0.1, 1.1, 1.0));
    camera.GetParams(3) = -1.01;
    EXPECT_TRUE(camera.IsBogusParams(0.1, 1.1, 1.0));
}

TEST(CImage, Default) 
{
    CImage image;
    EXPECT_EQ(image.GetImageID(), numeric_limits<size_t>::max());
    EXPECT_EQ(image.GetImageName(), "");
    EXPECT_EQ(image.GetCameraID(), numeric_limits<size_t>::max());
    EXPECT_FALSE(image.HasCamera());
    EXPECT_FALSE(image.IsRegistered());
    EXPECT_EQ(image.GetNumPoints2D(), 0);
    EXPECT_EQ(image.GetNumPoints3D(), 0);
    EXPECT_EQ(image.GetNumObservations(), 0);
    EXPECT_EQ(image.GetNumCorrespondences(), 0);
    EXPECT_EQ(image.GetNumVisiblePoints3D(), 0);
    EXPECT_EQ(image.GetPoint3DVisibilityScore(), 0);
    EXPECT_EQ(image.GetWorldToCamera().rotation.coeffs(), Eigen::Quaterniond::Identity().coeffs());
    EXPECT_EQ(image.GetWorldToCamera().translation, Eigen::Vector3d::Zero());
    auto a = image.GetWorldToCameraPrior().rotation.coeffs().array().isNaN();
    EXPECT_TRUE(image.GetWorldToCameraPrior().rotation.coeffs().array().isNaN().all());
    const CKeypoints points2D = image.GetKeypoints();
    EXPECT_EQ(points2D.size(), 0);
}

TEST(CImage, GetImageID)
{
    CImage image;
    EXPECT_EQ(image.GetImageID(), numeric_limits<size_t>::max());
    image.SetImageID(1);
    EXPECT_EQ(image.GetImageID(), 1);
}

TEST(CImage, Name)
{
    CImage image;
    EXPECT_EQ(image.GetImageName(), "");
    image.SetImageName("test1");
    EXPECT_EQ(image.GetImageName(), "test1");
    image.GetImageName() = "test2";
    EXPECT_EQ(image.GetImageName(), "test2");
}

TEST(CImage, CameraId)
{
    CImage image;
    EXPECT_EQ(image.GetCameraID(), numeric_limits<size_t>::max());
    image.SetCameraID(1);
    EXPECT_EQ(image.GetCameraID(), 1);
}

TEST(CImage, Registered)
{
    CImage image;
    EXPECT_FALSE(image.IsRegistered());
    image.SetRegistered(true);
    EXPECT_TRUE(image.IsRegistered());
    image.SetRegistered(false);
    EXPECT_FALSE(image.IsRegistered());
}

//TEST(CImage, NumPoints2D)
//{
//    CImage image;
//    EXPECT_EQ(image.GetNumPoints2D(), 0);
//    image.SetAllPoints2D(std::vector<Eigen::Vector2d>(10));
//    EXPECT_EQ(image.GetNumPoints2D(), 10);
//}
//
//TEST(CImage, NumPoints3D)
//{
//    CImage image;
//    image.SetAllPoints2D(std::vector<Eigen::Vector2d>(10));
//    EXPECT_EQ(image.GetNumPoints3D(), 0);
//    image.SetPoint3DForPoint2D(0, 0);
//    EXPECT_EQ(image.GetNumPoints3D(), 1);
//    image.SetPoint3DForPoint2D(0, 1);
//    image.SetPoint3DForPoint2D(1, 2);
//    EXPECT_EQ(image.GetNumPoints3D(), 2);
//}
//
//TEST(CImage, NumObservations)
//{
//    CImage image;
//    EXPECT_EQ(image.GetNumObservations(), 0);
//    image.SetNumObservations(10);
//    EXPECT_EQ(image.GetNumObservations(), 10);
//}
//
//TEST(CImage, NumCorrespondences)
//{
//    CImage image;
//    EXPECT_EQ(image.GetNumCorrespondences(), 0);
//    image.SetNumCorrespondences(10);
//    EXPECT_EQ(image.GetNumCorrespondences(), 10);
//}
//
//TEST(CImage, NumVisiblePoints3D)
//{
//    CImage image;
//    image.SetAllPoints2D(std::vector<Eigen::Vector2d>(10));
//    image.SetNumObservations(10);
//    CCamera camera;
//    camera.SetWidth(10);
//    camera.SetHeight(10);
//    image.Setup(camera);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 0);
//    image.IncrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 1);
//    image.IncrementCorrespondenceHasPoint3D(0);
//    image.IncrementCorrespondenceHasPoint3D(1);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 2);
//    image.DecrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 2);
//    image.DecrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 1);
//    image.DecrementCorrespondenceHasPoint3D(1);
//    EXPECT_EQ(image.GetNumVisiblePoints3D(), 0);
//}
//
//TEST(CImage, GetPoint3DVisibilityScore)
//{
//    CImage image;
//    std::vector<Eigen::Vector2d> points2D;
//    for (size_t i = 0; i < 4; ++i)
//    {
//        for (size_t j = 0; j < 4; ++j)
//        {
//            points2D.emplace_back(i, j);
//        }
//    }
//    image.SetAllPoints2D(points2D);
//    image.SetNumObservations(16);
//    CCamera camera;
//    camera.SetWidth(4);
//    camera.SetHeight(4);
//    image.Setup(camera);
//    Eigen::Matrix<size_t, Eigen::Dynamic, 1> scores(6, 1);
//    for (int i = 1; i <= 6; ++i)
//    {
//        scores(i - 1) = (1 << i) * (1 << i);
//    }
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), 0);
//    image.IncrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum());
//    image.IncrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum());
//    image.IncrementCorrespondenceHasPoint3D(1);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum() + scores.bottomRows(scores.size() - 1).sum());
//    image.IncrementCorrespondenceHasPoint3D(1);
//    image.IncrementCorrespondenceHasPoint3D(1);
//    image.IncrementCorrespondenceHasPoint3D(4);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum() + 2 * scores.bottomRows(scores.size() - 1).sum());
//    image.IncrementCorrespondenceHasPoint3D(4);
//    image.IncrementCorrespondenceHasPoint3D(5);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum() + 3 * scores.bottomRows(scores.size() - 1).sum());
//    image.DecrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum() + 3 * scores.bottomRows(scores.size() - 1).sum());
//    image.DecrementCorrespondenceHasPoint3D(0);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), scores.sum() + 2 * scores.bottomRows(scores.size() - 1).sum());
//    image.IncrementCorrespondenceHasPoint3D(2);
//    EXPECT_EQ(image.GetPoint3DVisibilityScore(), 2 * scores.sum() + 2 * scores.bottomRows(scores.size() - 1).sum());
//}
//
//TEST(CImage, Points2D) 
//{
//    CImage image;
//    std::vector<Eigen::Vector2d> points2D(10);
//    points2D[0] = Eigen::Vector2d(1.0, 2.0);
//    image.SetAllPoints2D(points2D);
//    EXPECT_EQ(image.GetAllPoints2D().size(), 10);
//    EXPECT_EQ(image.GetPoint2D(0).x(), 1.0);
//    EXPECT_EQ(image.GetPoint2D(0).y(), 2.0);
//    EXPECT_EQ(image.GetNumPoints3D(), 0);
//}
//
//TEST(CImage, Points2DWith3D)
//{
//    CImage image;
//    EXPECT_EQ(image.GetAllPoints2D().size(), 0);
//    std::vector<CPoint2D> points2D(10);
//    points2D[0].x() = 1.0;
//    points2D[0].y() = 2.0;
//    points2D[0].point3DID = 1;
//    image.SetAllPoints2D(points2D);
//    EXPECT_EQ(image.GetAllPoints2D().size(), 10);
//    EXPECT_EQ(image.GetPoint2D(0).x(), 1.0);
//    EXPECT_EQ(image.GetPoint2D(0).y(), 2.0);
//    EXPECT_EQ(image.GetNumPoints3D(), 1);
//}
//
//TEST(CImage, Points3D)
//{
//    CImage image;
//    image.SetAllPoints2D(std::vector<Eigen::Vector2d>(2));
//    EXPECT_FALSE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_FALSE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 0);
//    image.SetPoint3DForPoint2D(0, 0);
//    EXPECT_TRUE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_FALSE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 1);
//    EXPECT_TRUE(image.HasPoint3D(0));
//    image.SetPoint3DForPoint2D(0, 1);
//    EXPECT_TRUE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_FALSE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 1);
//    EXPECT_FALSE(image.HasPoint3D(0));
//    EXPECT_TRUE(image.HasPoint3D(1));
//    image.SetPoint3DForPoint2D(1, 0);
//    EXPECT_TRUE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_TRUE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 2);
//    EXPECT_TRUE(image.HasPoint3D(0));
//    EXPECT_TRUE(image.HasPoint3D(1));
//    image.ResetPoint3DForPoint2D(0);
//    EXPECT_FALSE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_TRUE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 1);
//    EXPECT_TRUE(image.HasPoint3D(0));
//    EXPECT_FALSE(image.HasPoint3D(1));
//    image.ResetPoint3DForPoint2D(1);
//    EXPECT_FALSE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_FALSE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 0);
//    EXPECT_FALSE(image.HasPoint3D(0));
//    EXPECT_FALSE(image.HasPoint3D(1));
//    image.ResetPoint3DForPoint2D(0);
//    EXPECT_FALSE(image.GetPoint2D(0).HasPoint3D());
//    EXPECT_FALSE(image.GetPoint2D(1).HasPoint3D());
//    EXPECT_EQ(image.GetNumPoints3D(), 0);
//    EXPECT_FALSE(image.HasPoint3D(0));
//    EXPECT_FALSE(image.HasPoint3D(1));
//}

TEST(CImage, ProjectionCenter)
{
    CImage image;
    EXPECT_EQ(image.GetProjectionCenter(), Eigen::Vector3d::Zero());
}

TEST(CImage, ViewingDirection)
{
    CImage image;
    EXPECT_EQ(image.GetViewDirection(), Eigen::Vector3d(0, 0, 1));
}

TEST(Point2D, Default)
{
    CPoint2D point2D;
    EXPECT_EQ(point2D.GetXY(), Eigen::Vector2d::Zero());
    EXPECT_EQ(point2D.point3DID, numeric_limits<size_t>::max());
    EXPECT_FALSE(point2D.HasPoint3D());
}

TEST(Point2D, Point3DId)
{
    CPoint2D point2D;
    EXPECT_EQ(point2D.point3DID, numeric_limits<size_t>::max());
    EXPECT_FALSE(point2D.HasPoint3D());
    point2D.point3DID = 1;
    EXPECT_EQ(point2D.point3DID, 1);
    EXPECT_TRUE(point2D.HasPoint3D());
    point2D.point3DID = numeric_limits<size_t>::max();
    EXPECT_EQ(point2D.point3DID, numeric_limits<size_t>::max());
    EXPECT_FALSE(point2D.HasPoint3D());
}

TEST(Point3D, Default)
{
    CPoint3D point3D;
    EXPECT_EQ(point3D.GetX(), 0);
    EXPECT_EQ(point3D.GetY(), 0);
    EXPECT_EQ(point3D.GetZ(), 0);
    EXPECT_EQ(point3D.GetXYZ()[0], point3D.GetX());
    EXPECT_EQ(point3D.GetXYZ()[1], point3D.GetY());
    EXPECT_EQ(point3D.GetXYZ()[2], point3D.GetZ());
    EXPECT_EQ(point3D.GetColor()[0], 0);
    EXPECT_EQ(point3D.GetColor()[1], 0);
    EXPECT_EQ(point3D.GetColor()[2], 0);
    EXPECT_EQ(point3D.GetError(), -1.0);
    EXPECT_FALSE(point3D.HasError());
    EXPECT_EQ(point3D.GetTrack().GetTrackLength(), 0);
}

TEST(Point3D, XYZ)
{
    CPoint3D point3D;
    EXPECT_EQ(point3D.GetX(), 0);
    EXPECT_EQ(point3D.GetY(), 0);
    EXPECT_EQ(point3D.GetZ(), 0);
    EXPECT_EQ(point3D.GetXYZ()[0], point3D.GetX());
    EXPECT_EQ(point3D.GetXYZ()[1], point3D.GetY());
    EXPECT_EQ(point3D.GetXYZ()[2], point3D.GetZ());
    point3D.SetXYZ(Eigen::Vector3d(0.1, 0.2, 0.3));
    EXPECT_EQ(point3D.GetX(), 0.1);
    EXPECT_EQ(point3D.GetY(), 0.2);
    EXPECT_EQ(point3D.GetZ(), 0.3);
    EXPECT_EQ(point3D.GetXYZ()[0], point3D.GetX());
    EXPECT_EQ(point3D.GetXYZ()[1], point3D.GetY());
    EXPECT_EQ(point3D.GetXYZ()[2], point3D.GetZ());
    point3D.SetXYZ(Eigen::Vector3d(0.2, 0.3, 0.4));
    EXPECT_EQ(point3D.GetX(), 0.2);
    EXPECT_EQ(point3D.GetY(), 0.3);
    EXPECT_EQ(point3D.GetZ(), 0.4);
    EXPECT_EQ(point3D.GetXYZ()[0], point3D.GetX());
    EXPECT_EQ(point3D.GetXYZ()[1], point3D.GetY());
    EXPECT_EQ(point3D.GetXYZ()[2], point3D.GetZ());
}

TEST(Point3D, RGB)
{
    CPoint3D point3D;
    EXPECT_EQ(point3D.GetColor()[0], 0);
    EXPECT_EQ(point3D.GetColor()[1], 0);
    EXPECT_EQ(point3D.GetColor()[2], 0);
    point3D.SetColor(Eigen::Vector3ub(1, 2, 3));
    EXPECT_EQ(point3D.GetColor()[0], 1);
    EXPECT_EQ(point3D.GetColor()[1], 2);
    EXPECT_EQ(point3D.GetColor()[2], 3);
}

TEST(Point3D, Error)
{
    CPoint3D point3D;
    EXPECT_EQ(point3D.GetError(), -1.0);
    EXPECT_FALSE(point3D.HasError());
    point3D.SetError(1.0);
    EXPECT_EQ(point3D.GetError(), 1.0);
    EXPECT_TRUE(point3D.HasError());
    point3D.SetError(-1.0);
    EXPECT_EQ(point3D.GetError(), -1.0);
    EXPECT_FALSE(point3D.HasError());
}

TEST(Point3D, Track)
{
    CPoint3D point3D;
    EXPECT_EQ(point3D.GetTrack().GetTrackLength(), 0);
    point3D.SetTrack(CTrack());
    EXPECT_EQ(point3D.GetTrack().GetTrackLength(), 0);
    CTrack track;
    track.AddElement(0, 1);
    track.AddElement(0, 2);
    point3D.SetTrack(track);
    EXPECT_EQ(point3D.GetTrack().GetTrackLength(), 2);
    track.AddElement(0, 3);
    EXPECT_EQ(point3D.GetTrack().GetTrackLength(), 2);
}


TEST(HasPointPositiveDepth, Nominal)
{
    const CRigid3D cam_from_world(Eigen::Quaterniond::Identity(),
        Eigen::Vector3d::Zero());
    const Eigen::Matrix3x4d cam_from_world_mat = cam_from_world.ToMatrix();

    // In the image plane
    const bool check1 = HasPointPositiveDepth(cam_from_world_mat, Eigen::Vector3d(0, 0, 0));
    EXPECT_FALSE(check1);
    const bool check2 =
        HasPointPositiveDepth(cam_from_world_mat, Eigen::Vector3d(0, 2, 0));
    EXPECT_FALSE(check2);

    // Infront of camera
    const bool check3 =
        HasPointPositiveDepth(cam_from_world_mat, Eigen::Vector3d(0, 0, 1));
    EXPECT_TRUE(check3);

    // Behind camera
    const bool check4 =
        HasPointPositiveDepth(cam_from_world_mat, Eigen::Vector3d(0, 0, -1));
    EXPECT_FALSE(check4);
}

TEST(CTrackElement, Empty) 
{
    CTrackElement track_el;
    EXPECT_EQ(track_el.imageID, numeric_limits<size_t>::max());
    EXPECT_EQ(track_el.point2DIndex, numeric_limits<size_t>::max());
}

TEST(Track, Default) {
    CTrack track;
    EXPECT_EQ(track.GetTrackLength(), 0);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
}

TEST(Track, SetElements) {
    CTrack track;
    std::vector<CTrackElement> elements;
    elements.emplace_back(0, 1);
    elements.emplace_back(0, 2);
    track.SetAllElements(elements);
    EXPECT_EQ(track.GetTrackLength(), 2);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
    EXPECT_EQ(track.GetElement(0).imageID, 0);
    EXPECT_EQ(track.GetElement(0).point2DIndex, 1);
    EXPECT_EQ(track.GetElement(1).imageID, 0);
    EXPECT_EQ(track.GetElement(1).point2DIndex, 2);
    for (size_t i = 0; i < track.GetTrackLength(); ++i) {
        EXPECT_EQ(track.GetElement(i).imageID, track.GetAllElements()[i].imageID);
        EXPECT_EQ(track.GetElement(i).point2DIndex, track.GetAllElements()[i].point2DIndex);
    }
}

TEST(Track, AddElement) {
    CTrack track;
    track.AddElement(0, 1);
    track.AddElement(CTrackElement(0, 2));
    std::vector<CTrackElement> elements;
    elements.emplace_back(0, 1);
    elements.emplace_back(0, 2);
    track.AddElements(elements);
    EXPECT_EQ(track.GetTrackLength(), 4);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
    EXPECT_EQ(track.GetElement(0).imageID, 0);
    EXPECT_EQ(track.GetElement(0).point2DIndex, 1);
    EXPECT_EQ(track.GetElement(1).imageID, 0);
    EXPECT_EQ(track.GetElement(1).point2DIndex, 2);
    EXPECT_EQ(track.GetElement(2).imageID, 0);
    EXPECT_EQ(track.GetElement(2).point2DIndex, 1);
    EXPECT_EQ(track.GetElement(3).imageID, 0);
    EXPECT_EQ(track.GetElement(3).point2DIndex, 2);
    for (size_t i = 0; i < track.GetTrackLength(); ++i) {
        EXPECT_EQ(track.GetElement(i).imageID, track.GetAllElements()[i].imageID);
        EXPECT_EQ(track.GetElement(i).point2DIndex, track.GetAllElements()[i].point2DIndex);
    }
}

TEST(Track, DeleteElement) {
    CTrack track;
    track.AddElement(0, 1);
    track.AddElement(0, 2);
    track.AddElement(0, 3);
    track.AddElement(0, 3);
    EXPECT_EQ(track.GetTrackLength(), 4);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
    track.DeleteElement(0);
    EXPECT_EQ(track.GetTrackLength(), 3);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
    EXPECT_EQ(track.GetElement(0).imageID, 0);
    EXPECT_EQ(track.GetElement(0).point2DIndex, 2);
    EXPECT_EQ(track.GetElement(1).imageID, 0);
    EXPECT_EQ(track.GetElement(1).point2DIndex, 3);
    EXPECT_EQ(track.GetElement(2).imageID, 0);
    EXPECT_EQ(track.GetElement(2).point2DIndex, 3);
    track.DeleteElement(0, 3);
    EXPECT_EQ(track.GetTrackLength(), 1);
    EXPECT_EQ(track.GetAllElements().size(), track.GetTrackLength());
    EXPECT_EQ(track.GetElement(0).imageID, 0);
    EXPECT_EQ(track.GetElement(0).point2DIndex, 2);
}


TEST(Track, Compress) {
    CTrack track;
    track.AddElement(0, 1);
    track.AddElement(0, 2);
    track.AddElement(0, 3);
    track.AddElement(0, 3);
    EXPECT_EQ(track.GetAllElements().capacity(), 4);
    track.DeleteElement(0);
    track.DeleteElement(0);
    track.Compress();
    EXPECT_EQ(track.GetAllElements().capacity(), 2);
}

TEST(CTwoViewGeometry, Default) {
    CTwoViewGeometry two_view_geometry;
    EXPECT_EQ(two_view_geometry.type, CTwoViewGeometryType::CUndefined);
    EXPECT_EQ(two_view_geometry.F, Eigen::Matrix3d::Zero());
    EXPECT_EQ(two_view_geometry.E, Eigen::Matrix3d::Zero());
    EXPECT_EQ(two_view_geometry.H, Eigen::Matrix3d::Zero());
    EXPECT_EQ(two_view_geometry.image1ToImage2.rotation.coeffs(),
        Eigen::Quaterniond::Identity().coeffs());
    EXPECT_EQ(two_view_geometry.image1ToImage2.translation,
        Eigen::Vector3d::Zero());
    EXPECT_TRUE(two_view_geometry.inlierMatches.empty());
}

TEST(CTwoViewGeometry, Invert) {
    CTwoViewGeometry two_view_geometry;
    two_view_geometry.type = CTwoViewGeometryType::CCalibrated;
    two_view_geometry.F = two_view_geometry.E = two_view_geometry.H =
        Eigen::Matrix3d::Identity();
    two_view_geometry.image1ToImage2 =
        CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 1, 2));
    two_view_geometry.inlierMatches.resize(2);
    two_view_geometry.inlierMatches[0] = CSIFTMatch(0, 1);
    two_view_geometry.inlierMatches[1] = CSIFTMatch(2, 3);

    two_view_geometry.Invert();
    EXPECT_EQ(two_view_geometry.type, CTwoViewGeometryType::CCalibrated);
    EXPECT_TRUE(two_view_geometry.F.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.E.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.H.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.image1ToImage2.rotation.isApprox(
        Eigen::Quaterniond::Identity()));
    EXPECT_TRUE(two_view_geometry.image1ToImage2.translation.isApprox(
        Eigen::Vector3d(-0, -1, -2)));
    EXPECT_EQ(two_view_geometry.inlierMatches[0].point2DIndex1, 1);
    EXPECT_EQ(two_view_geometry.inlierMatches[0].point2DIndex2, 0);
    EXPECT_EQ(two_view_geometry.inlierMatches[1].point2DIndex1, 3);
    EXPECT_EQ(two_view_geometry.inlierMatches[1].point2DIndex2, 2);

    two_view_geometry.Invert();
    EXPECT_EQ(two_view_geometry.type, CTwoViewGeometryType::CCalibrated);
    EXPECT_TRUE(two_view_geometry.F.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.E.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.H.isApprox(Eigen::Matrix3d::Identity()));
    EXPECT_TRUE(two_view_geometry.image1ToImage2.rotation.isApprox(
        Eigen::Quaterniond::Identity()));
    EXPECT_TRUE(two_view_geometry.image1ToImage2.translation.isApprox(
        Eigen::Vector3d(0, 1, 2)));
    EXPECT_EQ(two_view_geometry.inlierMatches[0].point2DIndex1, 0);
    EXPECT_EQ(two_view_geometry.inlierMatches[0].point2DIndex2, 1);
    EXPECT_EQ(two_view_geometry.inlierMatches[1].point2DIndex1, 2);
    EXPECT_EQ(two_view_geometry.inlierMatches[1].point2DIndex2, 3);
}

TEST(CVisibilityPyramid, Default) {
    CVisibilityPyramid pyramid;
    EXPECT_EQ(pyramid.NumLevels(), 0);
    EXPECT_EQ(pyramid.Width(), 0);
    EXPECT_EQ(pyramid.Height(), 0);
    EXPECT_EQ(pyramid.Score(), 0);
}

TEST(CVisibilityPyramid, Score) {
    for (int num_levels = 1; num_levels < 8; ++num_levels) {
        Eigen::VectorXi scores(num_levels);
        size_t max_score = 0;
        for (int i = 1; i <= num_levels; ++i) {
            scores(i - 1) = (1 << i) * (1 << i);
            max_score += scores(i - 1) * scores(i - 1);
        }

        CVisibilityPyramid pyramid(static_cast<size_t>(num_levels), 4, 4);
        EXPECT_EQ(pyramid.NumLevels(), num_levels);
        EXPECT_EQ(pyramid.Width(), 4);
        EXPECT_EQ(pyramid.Height(), 4);
        EXPECT_EQ(pyramid.Score(), 0);
        EXPECT_EQ(pyramid.MaxScore(), max_score);

        EXPECT_EQ(pyramid.Score(), 0);
        pyramid.SetPoint(0, 0);
        EXPECT_EQ(pyramid.Score(), scores.sum());
        pyramid.SetPoint(0, 0);
        EXPECT_EQ(pyramid.Score(), scores.sum());
        pyramid.SetPoint(0, 1);
        EXPECT_EQ(pyramid.Score(),
            scores.sum() + scores.tail(scores.size() - 1).sum());
        pyramid.SetPoint(0, 1);
        pyramid.SetPoint(0, 1);
        pyramid.SetPoint(1, 0);
        EXPECT_EQ(pyramid.Score(),
            scores.sum() + 2 * scores.tail(scores.size() - 1).sum());
        pyramid.SetPoint(1, 0);
        pyramid.SetPoint(1, 1);
        EXPECT_EQ(pyramid.Score(),
            scores.sum() + 3 * scores.tail(scores.size() - 1).sum());
        pyramid.ResetPoint(0, 0);
        EXPECT_EQ(pyramid.Score(),
            scores.sum() + 3 * scores.tail(scores.size() - 1).sum());
        pyramid.ResetPoint(0, 0);
        EXPECT_EQ(pyramid.Score(),
            scores.sum() + 2 * scores.tail(scores.size() - 1).sum());
        pyramid.SetPoint(0, 2);
        EXPECT_EQ(pyramid.Score(),
            2 * scores.sum() + 2 * scores.tail(scores.size() - 1).sum());
    }
}


int main(int argc, char* argv[])
{
    CSIFTMatch match(1, 2);
    string a = match.WriteToString();


    CSIFTMatch match2;
    match2.ReadFromString(a);
    string b = match2.WriteToString();
    bool flag = a == b;


    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return 0;
}