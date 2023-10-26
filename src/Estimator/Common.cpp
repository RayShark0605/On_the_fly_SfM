#include "Common.h"

using namespace std;

void CenterAndNormalizeImagePoints(const vector<Eigen::Vector2d>& points, vector<Eigen::Vector2d>& normedPoints, Eigen::Matrix3d& originToResult)
{
    Check(!points.empty());

    const size_t numPoints = points.size();

    // Step 1. ��������
    Eigen::Vector2d centroid(0, 0);
    for (const Eigen::Vector2d& point : points) 
    {
        centroid += point;
    }
    centroid /= numPoints;

    // Step 2. ����RMS(Root mean square)����. �������е㵽���ĵ�ƽ������
    double RMS_MeanDist = 0;
    for (const Eigen::Vector2d& point : points) 
    {
        RMS_MeanDist += (point - centroid).squaredNorm();
    }
    RMS_MeanDist = sqrt(RMS_MeanDist / numPoints);

    // Step 3. �����һ������. ����ƽ�ƺ�������Ϣ
    const double normFactor = sqrt(2.0) / RMS_MeanDist;
    originToResult << normFactor, 0, -normFactor * centroid(0), 0, normFactor, -normFactor * centroid(1), 0, 0, 1;

    // Step 4. Ӧ�ù�һ������. �ù�����Ĺ�һ������ȥת����ǰ2D��, �õ�ת��֮��ĵ�
    normedPoints.resize(numPoints);
    for (size_t i = 0; i < numPoints; i++)
    {
        normedPoints[i] = (originToResult * points[i].homogeneous()).hnormalized();
    }
}
void ComputeSquaredSampsonError(const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, vector<double>& residuals)
{
    Check(points1.size() == points2.size());
    
    // Step 1. ��ʼ��
    const size_t numPoints = points1.size(); // ��Ӧ�������
    residuals.resize(numPoints); // ÿһ�Զ�Ӧ���Sampson����ƽ��

    // Step 2. ѭ������
    for (size_t i = 0; i < numPoints; i++)
    {
        const Eigen::Vector3d epipolarLine1 = E * points1[i].homogeneous(); // ��һ�����Ӧ�ļ���
        const Eigen::Vector3d point2Homogeneous = points2[i].homogeneous(); // �ڶ������Ӧ���������
        const double num = point2Homogeneous.dot(epipolarLine1); // Sampson������ʽ�еķ���(p2^T * E * p1)
        const Eigen::Vector4d denom(point2Homogeneous.dot(E.col(0)), point2Homogeneous.dot(E.col(1)), epipolarLine1.x(), epipolarLine1.y());
        residuals[i] = num * num / denom.squaredNorm();
    }
}
void ComputeSquaredReprojectionError(const vector<Eigen::Vector2d>& points2D, const vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& worldToCamera, vector<double>& residuals)
{
    Check(points2D.size() == points3D.size());

    // Step 1. ��ʼ��
    const size_t numPoints = points2D.size();
    residuals.resize(numPoints);

    // �������е�Լ������
    for (size_t i = 0; i < numPoints; i++)
    {
        const Eigen::Vector3d point3DInCamera = worldToCamera * points3D[i].homogeneous(); // ��3D�����������ϵת�����������ϵ
        if (point3DInCamera.z() > numeric_limits<double>::epsilon()) // ���3D�����������ϵ�����zֵ����0, ��ô����Ϊ���������ǰ��
        {
            residuals[i] = (point3DInCamera.hnormalized() - points2D[i]).squaredNorm(); // ����õ�Ĺ�һ���������Ӧ2D��֮��ľ���ƽ��
        }
        else
        {
            residuals[i] = numeric_limits<double>::max();
        }
    }
}
