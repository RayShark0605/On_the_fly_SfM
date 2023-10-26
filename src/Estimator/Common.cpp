#include "Common.h"

using namespace std;

void CenterAndNormalizeImagePoints(const vector<Eigen::Vector2d>& points, vector<Eigen::Vector2d>& normedPoints, Eigen::Matrix3d& originToResult)
{
    Check(!points.empty());

    const size_t numPoints = points.size();

    // Step 1. 计算质心
    Eigen::Vector2d centroid(0, 0);
    for (const Eigen::Vector2d& point : points) 
    {
        centroid += point;
    }
    centroid /= numPoints;

    // Step 2. 计算RMS(Root mean square)距离. 计算所有点到质心的平均距离
    double RMS_MeanDist = 0;
    for (const Eigen::Vector2d& point : points) 
    {
        RMS_MeanDist += (point - centroid).squaredNorm();
    }
    RMS_MeanDist = sqrt(RMS_MeanDist / numPoints);

    // Step 3. 构造归一化矩阵. 包含平移和缩放信息
    const double normFactor = sqrt(2.0) / RMS_MeanDist;
    originToResult << normFactor, 0, -normFactor * centroid(0), 0, normFactor, -normFactor * centroid(1), 0, 0, 1;

    // Step 4. 应用归一化矩阵. 用构造出的归一化矩阵去转换当前2D点, 得到转换之后的点
    normedPoints.resize(numPoints);
    for (size_t i = 0; i < numPoints; i++)
    {
        normedPoints[i] = (originToResult * points[i].homogeneous()).hnormalized();
    }
}
void ComputeSquaredSampsonError(const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, vector<double>& residuals)
{
    Check(points1.size() == points2.size());
    
    // Step 1. 初始化
    const size_t numPoints = points1.size(); // 对应点的数量
    residuals.resize(numPoints); // 每一对对应点的Sampson误差的平方

    // Step 2. 循环计算
    for (size_t i = 0; i < numPoints; i++)
    {
        const Eigen::Vector3d epipolarLine1 = E * points1[i].homogeneous(); // 第一个点对应的极线
        const Eigen::Vector3d point2Homogeneous = points2[i].homogeneous(); // 第二个点对应的齐次坐标
        const double num = point2Homogeneous.dot(epipolarLine1); // Sampson误差计算式中的分子(p2^T * E * p1)
        const Eigen::Vector4d denom(point2Homogeneous.dot(E.col(0)), point2Homogeneous.dot(E.col(1)), epipolarLine1.x(), epipolarLine1.y());
        residuals[i] = num * num / denom.squaredNorm();
    }
}
void ComputeSquaredReprojectionError(const vector<Eigen::Vector2d>& points2D, const vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& worldToCamera, vector<double>& residuals)
{
    Check(points2D.size() == points3D.size());

    // Step 1. 初始化
    const size_t numPoints = points2D.size();
    residuals.resize(numPoints);

    // 遍历所有点对计算误差
    for (size_t i = 0; i < numPoints; i++)
    {
        const Eigen::Vector3d point3DInCamera = worldToCamera * points3D[i].homogeneous(); // 将3D点从世界坐标系转换到相机坐标系
        if (point3DInCamera.z() > numeric_limits<double>::epsilon()) // 如果3D点在相机坐标系坐标的z值大于0, 那么就认为它在相机的前面
        {
            residuals[i] = (point3DInCamera.hnormalized() - points2D[i]).squaredNorm(); // 计算该点的归一化坐标与对应2D点之间的距离平方
        }
        else
        {
            residuals[i] = numeric_limits<double>::max();
        }
    }
}
