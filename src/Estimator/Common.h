#pragma once
#define EIGEN_USE_MKL_ALL
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Sampler.h"
#include "SupportMeasurer.h"
#include "../Base/Base.h"


// ���Ļ��͹�һ��ͼ���ϵĵ�. 
// �ᾭ�������任: 1. ���Ļ�: �µ�����ϵͳ��ԭ�㽫λ��ͼ��������. 2. ��һ��: ʹ�õ㵽������ϵͳԭ���ƽ������Ϊsqrt(2)
//
// @param points            ������.
// @param normedPoints      ���: ת����ĵ�����.
// @param originToResult    ���: 3��3��ת������.
void CenterAndNormalizeImagePoints(const std::vector<Eigen::Vector2d>& points, std::vector<Eigen::Vector2d>& normedPoints, Eigen::Matrix3d& originToResult);

// ����һ���Ӧ��͸����Ļ���������ʾ���Ĳв�. �в����ΪSampson����ƽ��
//
// @param points1     ��һ���Ӧ��.
// @param points2     �ڶ����Ӧ��.
// @param E           3��3�Ļ���������ʾ���
// @param residuals   ���: �в�����.
void ComputeSquaredSampsonError(const std::vector<Eigen::Vector2d>& points1, const std::vector<Eigen::Vector2d>& points2, const Eigen::Matrix3d& E, std::vector<double>& residuals);

// ����һ��2Dͼ������֮��Ӧ��3D���ƽ����ͶӰ���. ���ĳ��3D��������ĺ���, ��ô����Ӧ��ƽ����ͶӰ��������Ϊstd::numeric_limits<double>::max()
//
// @param points2D      ��һ����2Dͼ���
// @param points3D      3D��
// @param proj_matrix   3��4��ͶӰ����
// @param residuals     ���: ƽ����ͶӰ�������
void ComputeSquaredReprojectionError(const std::vector<Eigen::Vector2d>& points2D, const std::vector<Eigen::Vector3d>& points3D, const Eigen::Matrix3x4d& worldToCamera, std::vector<double>& residuals);
