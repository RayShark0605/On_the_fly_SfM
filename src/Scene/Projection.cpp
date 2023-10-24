#include "Projection.h"

using namespace std;

double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera)
{
    const Eigen::Vector3d point3D_Camera = worldToCamera * point3D;

    // 检查3D点是否在相机前方
    if (point3D_Camera.z() < numeric_limits<double>::epsilon()) 
    {
        return numeric_limits<double>::max();
    }
    float imagePointX, imagePointY;
    const Eigen::Vector2d cameraPoint = point3D_Camera.hnormalized();
    camera.CameraToImage(cameraPoint(0), cameraPoint(1), imagePointX, imagePointY);
    return (Eigen::Vector2d(imagePointX, imagePointY) - point2D).squaredNorm();
}
double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera)
{
    const double projectZ = worldToCamera.row(2).dot(point3D.homogeneous());

    // 检查3D点是否在相机前方
    if (projectZ < numeric_limits<double>::epsilon())
    {
        return numeric_limits<double>::max();
    }

    const double projectX = worldToCamera.row(0).dot(point3D.homogeneous());
    const double projectY = worldToCamera.row(1).dot(point3D.homogeneous());
    const double invProjectZ = 1.0 / projectZ;

    Eigen::Vector2d CameraPoint(invProjectZ * projectX, invProjectZ * projectY);
    float imagePointX, imagePointY;
    camera.CameraToImage(CameraPoint(0), CameraPoint(1), imagePointX, imagePointY);
    const Eigen::Vector2d projectedPoint2D(imagePointX, imagePointY);
    return (projectedPoint2D - point2D).squaredNorm();
}

double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera, const CCamera& camera)
{
    float cameraPointX, cameraPointY;
    camera.ImageToCamera(point2D(0), point2D(1), cameraPointX, cameraPointY);
    return CalculateNormalizedAngularError(Eigen::Vector2d(cameraPointX, cameraPointY), point3D, worldToCamera);
}
double CalculateAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera, const CCamera& camera)
{
    float cameraPointX, cameraPointY;
    camera.ImageToCamera(point2D(0), point2D(1), cameraPointX, cameraPointY);
    return CalculateNormalizedAngularError(Eigen::Vector2d(cameraPointX, cameraPointY), point3D, worldToCamera);
}

double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const CRigid3D& worldToCamera)
{
    const Eigen::Vector3d ray1 = point2D.homogeneous();
    const Eigen::Vector3d ray2 = worldToCamera * point3D;
    return acos(ray1.normalized().transpose() * ray2.normalized());
}
double CalculateNormalizedAngularError(const Eigen::Vector2d& point2D, const Eigen::Vector3d& point3D, const Eigen::Matrix3x4d& worldToCamera)
{
    const Eigen::Vector3d ray1 = point2D.homogeneous();
    const Eigen::Vector3d ray2 = worldToCamera * point3D.homogeneous();
    return acos(ray1.normalized().transpose() * ray2.normalized());
}

bool HasPointPositiveDepth(const Eigen::Matrix3x4d& worldToCamera, const Eigen::Vector3d& point3D)
{
    return worldToCamera.row(2).dot(point3D.homogeneous()) >= numeric_limits<double>::epsilon();
}
