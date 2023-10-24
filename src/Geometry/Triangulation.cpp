#include "Triangulation.h"

using namespace std;

Eigen::Vector3d TriangulatePoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
{
	// Step 1. ����˫�Ӽ���֮������ǲ�����������������A
	Eigen::Matrix4d A;
	A.row(0) = point1(0) * worldToCamera1.row(2) - worldToCamera1.row(0);
	A.row(1) = point1(1) * worldToCamera1.row(2) - worldToCamera1.row(1);
	A.row(2) = point2(0) * worldToCamera2.row(2) - worldToCamera2.row(0);
	A.row(3) = point2(1) * worldToCamera2.row(2) - worldToCamera2.row(1);

	// Step 2. ��A������ֵ�ֽ�, ��������������������
	Eigen::JacobiSVD<Eigen::Matrix4d> svd(A, Eigen::ComputeFullV);

	// Step 3. ���һ��������������Ӧ�ھ���A����С����ֵ, �����ǲ���������, ��������ͱ�ʾ��ά��. ͨ��hnormalizedתΪ�������
	return svd.matrixV().col(3).hnormalized();
}
Eigen::Vector3d TriangulateOptimalPoint(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
{
	const CRigid3D worldToCamera1_Rigid3D(Eigen::Quaterniond(worldToCamera1.leftCols<3>()), worldToCamera1.col(3));
	const CRigid3D worldToCamera2_Rigid3D(Eigen::Quaterniond(worldToCamera2.leftCols<3>()), worldToCamera2.col(3));
	const CRigid3D camera1ToCamera2_Rigid3D = worldToCamera2_Rigid3D * worldToCamera1_Rigid3D.Inverse();
	const Eigen::Matrix3d E = PoseToEssentialMatrix(camera1ToCamera2_Rigid3D);

	Eigen::Vector2d optimalPoint1;
	Eigen::Vector2d optimalPoint2;
	FindOptimalImageObservations(E, point1, point2, optimalPoint1, optimalPoint2);
	return TriangulatePoint(worldToCamera1, worldToCamera2, optimalPoint1, optimalPoint2);
}
Eigen::Vector3d TriangulateMultiViewPoint(const vector<Eigen::Matrix3x4d>& worldToCameras, const vector<Eigen::Vector2d>& points2D)
{
	CHECK(worldToCameras.size() == points2D.size());

	// Step 1. ��ʼ��4��4�ľ���A, ���������ۻ�������ͼ����Ϣ
	Eigen::Matrix4d A = Eigen::Matrix4d::Zero();

	// Step 2. �������е���ͼ��ÿ����ͼ�ж�Ӧ��2D��
	for (size_t i = 0; i < points2D.size(); i++)
	{
		// Step 2.1. ��2D��ת��Ϊ�������, �����й�һ��
		const Eigen::Vector3d point = points2D[i].homogeneous().normalized();

		// Step 2.2. ���������(������): ԭʼͶӰ�����ȥ�������ͶӰ����
		const Eigen::Matrix3x4d term = worldToCameras[i] - point * point.transpose() * worldToCameras[i];

		// Step 2.3. ���������ƽ�����
		A += term.transpose() * term;
	}

	// Step 3. ʹ������ֵ�ֽ�������С��������. ����ʹ�õ����԰����(self-adjoint)������ֵ�����, ����ֵ���ȶ�
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigenSolver(A);

	// Step 4. ��������С����ֵ��Ӧ����������, ���������ŵ�3D��λ��
	return eigenSolver.eigenvectors().col(0).hnormalized();
}
vector<Eigen::Vector3d> TriangulatePoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2)
{
	CHECK(points1.size() == points2.size());
	vector<Eigen::Vector3d> points3D(points1.size());
	for (size_t i = 0; i < points3D.size(); i++) 
	{
		points3D[i] = TriangulatePoint(worldToCamera1, worldToCamera2, points1[i], points2[i]);
	}
	return points3D;
}
vector<Eigen::Vector3d> TriangulateOptimalPoints(const Eigen::Matrix3x4d& worldToCamera1, const Eigen::Matrix3x4d& worldToCamera2, const vector<Eigen::Vector2d>& points1, const vector<Eigen::Vector2d>& points2)
{
	CHECK(points1.size() == points2.size());
	vector<Eigen::Vector3d> points3D(points1.size());
	for (size_t i = 0; i < points3D.size(); i++)
	{
		points3D[i] = TriangulateOptimalPoint(worldToCamera1, worldToCamera2, points1[i], points2[i]);
	}
	return points3D;
}
double CalculateTriangulationAngle(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const Eigen::Vector3d& point3D)
{
	const double baselineLengthSquared = (projectionCenter1 - projectionCenter2).squaredNorm();

	const double rayLengthSquared1 = (point3D - projectionCenter1).squaredNorm();
	const double rayLengthSquared2 = (point3D - projectionCenter2).squaredNorm();

	// ʹ�����Ҷ�����㽻���
	const double denominator = 2.0 * sqrt(rayLengthSquared1 * rayLengthSquared2);
	if (denominator == 0.0)
	{
		return 0.0;
	}
	const double nominator = rayLengthSquared1 + rayLengthSquared2 - baselineLengthSquared;
	const double angle = abs(acos(nominator / denominator));

	// ���ǲ��������(�Ƚ�Զ��ͬ����)�Ͷ۽�(�ȽϽ���ͬ����)��������ǲ��ȶ���, ��������Ҫ���������ཻ�����γɵ����ֽǶ�(��ǺͶ۽�)�еĽ�Сֵ
	return min(angle, M_PI - angle);
}
vector<double> CalculateTriangulationAngles(const Eigen::Vector3d& projectionCenter1, const Eigen::Vector3d& projectionCenter2, const vector<Eigen::Vector3d>& points3D)
{
	const double baselineSquared = (projectionCenter1 - projectionCenter2).squaredNorm();
	vector<double> angles(points3D.size());
	for (size_t i = 0; i < points3D.size(); i++)
	{
		const double rayLengthSquared1 = (points3D[i] - projectionCenter1).squaredNorm();
		const double rayLengthSquared2 = (points3D[i] - projectionCenter2).squaredNorm();

		const double denominator = 2.0 * sqrt(rayLengthSquared1 * rayLengthSquared2);
		if (denominator == 0.0)
		{
			angles[i] = 0.0;
			continue;
		}
		const double nominator = rayLengthSquared1 + rayLengthSquared2 - baselineSquared;
		const double angle = abs(acos(nominator / denominator));
		angles[i] = min(angle, M_PI - angle);
	}
	return angles;
}
bool EstimateTriangulation(const CEstimateTriangulationOptions& options, const vector<CTriangulationPoint>& points, const vector<CTriangulationPose>& poses, vector<char>& inlierMask, Eigen::Vector3d& XYZ)
{
	CHECK(points.size() >= 2 && points.size() == poses.size());
	options.Check();

	CTriangulationEstimator triangulationEstimator(options.minTriAngle, options.residualType);
	CTriangulationEstimator triangulationEstimator_Local(options.minTriAngle, options.residualType);
	CInlierSupportMeasurer inlierSupportMeasurer;
	CCombinationSampler combinationSampler(triangulationEstimator.minNumSamples);
	CLORANSAC loransac(options.ransacOptions, &triangulationEstimator, &triangulationEstimator_Local, &inlierSupportMeasurer, &combinationSampler);
	const CRANSACReport report = loransac.Estimate<CTriangulationPoint, CTriangulationPose, Eigen::Vector3d>(points, poses);
	if (!report.isSuccess) 
	{
		return false;
	}
	inlierMask = report.inlierMask;
	XYZ = report.model;
	return true;
}
