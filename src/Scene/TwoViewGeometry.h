#pragma once
#define EIGEN_USE_MKL_ALL
#include <algorithm>
#include <vector>

#include <tbb/tbb.h>
#include <Eigen/Core>

#include "Point2D.h"
#include "Point3D.h"
#include "Camera.h"
#include "../Base/Options.h"
#include "../Base/Math.h"
#include "../Geometry/Rigid3D.h"
#include "../Geometry/Triangulation.h"
#include "../Geometry/EssentialMatrix.h"
#include "../Geometry/HomographyMatrix.h"

// ˫�Ӽ��ι�ϵ����
enum CTwoViewGeometryType :size_t
{
	CUndefined = 0,              // ˫�Ӽ��ι�ϵ��δ����
	CDegenerate = 1,             // �˻���˫�Ӽ��ι�ϵ, �޷��ƶϳ��Ƚ��ļ���ģ��. ԭ��: ˫��ͼ֮���ص����򲻹�, ƥ���(�ڵ�)����
	CCalibrated = 2,             // ��У׼��˫�Ӽ��ι�ϵ, ʹ�ñ��ʾ�������, ����ڲ���֪, �����ƶϳ�׼ȷ��˫�Ӽ��ι�ϵ
	CUncalibrated = 3,           // δУ׼��˫�Ӽ��ι�ϵ, ʹ�û�����������. ԭ��: ����ڲ�δ֪��׼ȷ
	CPlanar = 4,                 // ƽ��˫�Ӽ��ι�ϵ, ����ƽ�ƹ�ϵ, ʹ�õ�Ӧ����
	CPanoramic = 5,              // ��ƽ�ƹ�ϵ�Ĵ���ת˫�Ӽ��ι�ϵ, ʹ�õ�Ӧ����
	CPlanarOrPanoramic = 6,      // ƽ�����ת˫�Ӽ��ι�ϵ, ʹ�õ�Ӧ����
	CWatermark = 7,              // Ӱ���Ե�Ĵ���άƽ��. һ����ͼ���Ե��ˮӡ���߱�ǩ
	CMultiple = 8,               // ���϶���˫�Ӽ������͵Ļ��. ����: ��������Щ������CCalibrated����, Ҳ��һЩ����(����ǽ��)��CPlanar����
	CTwoViewGeometryTypeNum = 9, // ����9��˫�Ӽ�������
};

class CTwoViewGeometry final
{
public:
	size_t type = CTwoViewGeometryType::CUndefined; // ˫�Ӽ��ι�ϵ����
	Eigen::Matrix3d E = Eigen::Matrix3d::Zero();    // ���ʾ���
	Eigen::Matrix3d F = Eigen::Matrix3d::Zero();    // ��������
	Eigen::Matrix3d H = Eigen::Matrix3d::Zero();    // ��Ӧ����
	CRigid3D image1ToImage2;                        // Ӱ��1��Ӱ��2�����λ��(����任)
	CSIFTMatches inlierMatches;                     // �ڵ�
	double meanTriAngle = -1;                       // ���ǻ�����ǵ���λ��

	// ��ת˫�Ӽ��ι�ϵ��ƥ�佻����Ӱ���(����ǰ˫�Ӽ��ι�ϵ����Ϊ: ��Ӱ��2��Ӱ��1�����λ��). ���������ı��Ա����.
	void Invert();

	// ��Զ���. ����˫�Ӽ��ε����λ��
	bool EstimateRelativePose(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2);

	// ����������Ƿ�������齹��, ��У׼��δУ׼��Ӱ����й���˫�Ӽ���
	void Estimate(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, const CSIFTMatches& matches, const CTwoViewGeometryOptions& options);

	// ��У׼��Ӱ����й���˫�Ӽ���
	void EstimateCalibrated(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, const CSIFTMatches& matches, const CTwoViewGeometryOptions& options);

	// ��δУ׼��Ӱ����й���˫�Ӽ���
	void EstimateUncalibrated(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, const CSIFTMatches& matches, const CTwoViewGeometryOptions& options);

	// ��У׼��Ӱ����й��Ƶ�Ӧ����
	void EstimateCalibratedHomography(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, const CSIFTMatches& matches, const CTwoViewGeometryOptions& options);
	
	// ���ƶ�����˫�Ӽ���. ��ƥ�����Ƴ���ǰ���ڵ㼯���ݹ���ƶ������, ֱ���Ҳ����㹻���ڵ�Ϊֹ. ����ܹ��Ƴ��������, ��ô�ڵ�ƥ�佫����������, �������ͽ�ΪCMultiple
	// �����ھ��и��ӳ���(��������ƶ�����)�Լ�������Ӱ��
	void EstimateMultipleType(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, const CSIFTMatches& matches, const CTwoViewGeometryOptions& options);

	// ����ڵ�ƥ���Ƿ���ˮӡ����. ˮӡ��ͼ��ı�Ե��������ƽ��
	bool DetectWaterMark(const CCamera& camera1, const std::vector<Eigen::Vector2d>& points1, const CCamera& camera2, const std::vector<Eigen::Vector2d>& points2, size_t numInliers, const std::vector<char>& inlierMask, const CTwoViewGeometryOptions& options);
};