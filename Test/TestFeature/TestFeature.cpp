#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <opencv2/xfeatures2d.hpp>
#include <immintrin.h>
#include <Eigen/Dense>
#include <random>
#include <cuda_runtime.h>

#include "../../src/Feature/Common.h"
#include "../../src/Feature/FeatureExtractor.h"
#include "../../src/Feature/FeatureMatcher.h"
#include "../../src/Scene/TwoViewGeometry.h"

using namespace std;


Eigen::MatrixXf ComputeDistances(const Eigen::Matrix<unsigned char, Eigen::Dynamic, 128, Eigen::RowMajor>& a, const Eigen::Matrix<unsigned char, Eigen::Dynamic, 128, Eigen::RowMajor>& b)
{

    //OpenMP 500多ms
    int m = a.rows();
    int n = b.rows();

    Eigen::MatrixXf distances(m, n);

    _Pragma("omp parallel for")
    for (int i = 0; i < m;i++)
    {
        for (int j = 0; j < n; j++)
        {
            distances(i, j) = (a.row(i).cast<float>() - b.row(j).cast<float>()).norm();
        }
    }

    return distances;
    

}






int main(int argc, char* argv[])
{
    /*const int Dim = 128;
    const int Rows1 = 8000;
    const int Rows2 = 9000;
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Dim, Eigen::RowMajor> matrix1(Rows1, Dim);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Dim, Eigen::RowMajor> matrix2(Rows2, Dim);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 255);

    for (int i = 0; i < Rows1; i++)
        for (int j = 0; j < Dim; j++)
            matrix1(i, j) = static_cast<unsigned char>(dis(gen));

    for (int i = 0; i < Rows2; i++)
        for (int j = 0; j < Dim; j++)
            matrix2(i, j) = static_cast<unsigned char>(dis(gen));




    Eigen::MatrixXf distancesTrue(Rows1, Rows2);
    for (int i = 0; i < Rows1; i++)
    {
        for (int j = 0; j < Rows2; j++)
        {
            distancesTrue(i, j) = (matrix1.row(i).cast<float>() - matrix2.row(j).cast<float>()).norm();
        }
    }
    
    cout << "开始计算!" << endl;
    auto start = chrono::high_resolution_clock::now();
    Eigen::MatrixXf distancesCalculate = ComputeDistances(matrix1, matrix2);
    auto stop = chrono::high_resolution_clock::now();
    cout << "耗时: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << "ms" << endl;

    float tolerance = 1e-3;
    if (distancesCalculate.isApprox(distancesTrue, tolerance))
    {
        cout << "计算正确!" << endl;
    }
    else
    {
        cout << "计算错误!" << endl;
        cout << distancesCalculate << endl << endl << endl;
        cout << distancesTrue << endl << endl << endl;
        cout << matrix1 << endl;
    }*/


	int index = SetBestCudaDevice();
    CSIFTExtractionOptions SIFTExtractionOptions;
    CSIFTGPUExtractor SIFTGPUExtractor(SIFTExtractionOptions);

    string image1Path = "C:/Users/XiaRui/Desktop/PTRS_Project/images/DSC00013.JPG";
    CKeypoints keypoints1;
    CSIFTDescriptors descriptors1;
    CPartedTime p1(1);
    if (!SIFTGPUExtractor.Extract(image1Path, keypoints1, descriptors1))
    {
        cout << "提取失败!" << endl;
    }
    p1.~CPartedTime();


    CPartedTime p2(2);
    string image2Path = "C:/Users/XiaRui/Desktop/PTRS_Project/images/DSC00014.JPG";
    CKeypoints keypoints2;
    CSIFTDescriptors descriptors2;
    if (!SIFTGPUExtractor.Extract(image2Path, keypoints2, descriptors2))
    {
        cout << "提取失败!" << endl;
    }
    p2.~CPartedTime();

    CPartedTime p3(3);
    CSIFTMatchingOptions SIFTMatchingOptions;
    CSIFTCPUMatcher SIFTCPUMatcher(SIFTMatchingOptions);
    CSIFTMatches matches = SIFTCPUMatcher.Match(descriptors1, descriptors2);
    cout << "匹配数量: " << matches.size() << endl;
    p3.~CPartedTime();

    CPartedTime p4(4);
    CTwoViewGeometry twoViewGeometry;
    CCamera camera;
    camera.SetCameraModelName("Simple Radial");
    vector<Eigen::Vector2d> points1(keypoints1.size()), points2(keypoints2.size());
    for (size_t i = 0; i < keypoints1.size(); i++)
    {
        points1[i] = Eigen::Vector2d(keypoints1[i].pt.x, keypoints1[i].pt.y);
    }
    for (size_t i = 0; i < keypoints2.size(); i++)
    {
        points2[i] = Eigen::Vector2d(keypoints2[i].pt.x, keypoints2[i].pt.y);
    }
    CTwoViewGeometryOptions twoViewGeometryOptions;
    twoViewGeometryOptions.isComputeRelativePose = true;
    twoViewGeometry.Estimate(camera, points1, camera, points2, matches, twoViewGeometryOptions);
    cout << "内匹配数量: " << twoViewGeometry.inlierMatches.size() << endl;
    p4.~CPartedTime();

    cv::Mat img1 = cv::imread(image1Path, cv::IMREAD_COLOR);
    cv::Mat img2 = cv::imread(image2Path, cv::IMREAD_COLOR);
    vector<cv::DMatch> dmatches;
    for (size_t i = 0; i < matches.size(); i++)
    {
        dmatches.push_back(cv::DMatch(matches[i].point2DIndex1, matches[i].point2DIndex2, 0));
    }
    vector<cv::KeyPoint> keypointsCV1(keypoints1.size());
    for (size_t i = 0; i < keypoints1.size(); i++)
    {
        keypointsCV1[i] = keypoints1[i];
    }
    vector<cv::KeyPoint> keypointsCV2(keypoints2.size());
    for (size_t i = 0; i < keypoints2.size(); i++)
    {
        keypointsCV2[i] = keypoints2[i];
    }
    cv::Mat outImg;
    cv::Size size(img1.cols + img2.cols, std::max(img1.rows, img2.rows));
    outImg.create(size, CV_MAKETYPE(img1.depth(), 3));
    outImg = cv::Scalar::all(0);

    // Copy images to output
    img1.copyTo(outImg(cv::Rect(0, 0, img1.cols, img1.rows)));
    img2.copyTo(outImg(cv::Rect(img1.cols, 0, img2.cols, img2.rows)));

    // Draw matches
    for (const auto& match : twoViewGeometry.inlierMatches)
    {
        cv::Point2f pt1 = keypoints1[match.point2DIndex1].pt;
        cv::Point2f pt2 = keypoints2[match.point2DIndex2].pt;
        pt2.x += img1.cols;  // Shift x coordinate for second image

        cv::line(outImg, pt1, pt2, cv::Scalar(0, 255, 0), 4);
    }
    cv::imwrite("1111.bmp", outImg);




    /*cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches, outImg,);
    cv::imwrite("1111.bmp", outImg);*/
    /*vector<cv::DMatch> dmatches2;
    for (size_t i = 0; i < twoViewGeometry.inlierMatches.size(); i++)
    {
        dmatches2.push_back(cv::DMatch(twoViewGeometry.inlierMatches[i].point2DIndex1, twoViewGeometry.inlierMatches[i].point2DIndex2, 0));
    }
    cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches2, outImg);
    cv::imwrite("7twoViewGeometry.bmp", outImg);*/


    ///*CPartedTime p5(5);
    //SIFTCPUMatcher.MatchGuided(keypoints1, keypoints2, descriptors1, descriptors2, twoViewGeometry);
    //cout << "指导匹配后, 内匹配数量: " << twoViewGeometry.inlierMatches.size() << endl;
    //p5.~CPartedTime();
    //vector<cv::DMatch> dmatches3;
    //for (size_t i = 0; i < twoViewGeometry.inlierMatches.size(); i++)
    //{
    //    dmatches3.push_back(cv::DMatch(twoViewGeometry.inlierMatches[i].point2DIndex1, twoViewGeometry.inlierMatches[i].point2DIndex2, 0));
    //}
    //cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches3, outImg);
    //cv::imwrite("8guideMatch.bmp", outImg);

    //cout << endl << twoViewGeometry.image1ToImage2.ToMatrix() << endl;
    //cout << endl << twoViewGeometry.E << endl;
    //cout << endl << twoViewGeometry.F << endl;
    //cout << endl << twoViewGeometry.H << endl;
    //cout << endl << twoViewGeometry.meanTriAngle << endl;*/
	
    //CSIFTGPUMatcher SIFTGPUMatcher(SIFTMatchingOptions);
    //auto start = chrono::high_resolution_clock::now();
    //CSIFTMatches gpuMatches = SIFTGPUMatcher.Match(descriptors1, descriptors2);
    //auto end = chrono::high_resolution_clock::now();
    //cout << "GPU SIFT匹配时间: " << chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << endl;
    //cout << "GPU SIFT匹配数: " << gpuMatches.size() << endl;
    //vector<cv::DMatch> dmatches4;
    //for (size_t i = 0; i < gpuMatches.size(); i++)
    //{
    //    dmatches4.push_back(cv::DMatch(gpuMatches[i].point2DIndex1, gpuMatches[i].point2DIndex2, 0));
    //}
    //cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches4, outImg);
    //cv::imwrite("9gpuMatch.bmp", outImg);

    //start = chrono::high_resolution_clock::now();
    //CTwoViewGeometry gpuTwoViewGeometry;
    //gpuTwoViewGeometry.Estimate(camera, points1, camera, points2, gpuMatches, twoViewGeometryOptions);
    //end = chrono::high_resolution_clock::now();
    //cout << "GPU 双视几何估计时间: " << chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << endl;
    //cout << "GPU内匹配数量: " << gpuTwoViewGeometry.inlierMatches.size() << endl;
    //vector<cv::DMatch> dmatches5;
    //for (size_t i = 0; i < gpuTwoViewGeometry.inlierMatches.size(); i++)
    //{
    //    dmatches5.push_back(cv::DMatch(gpuTwoViewGeometry.inlierMatches[i].point2DIndex1, gpuTwoViewGeometry.inlierMatches[i].point2DIndex2, 0));
    //}
    //cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches5, outImg);
    //cv::imwrite("10GPUTwoViewGeometry.bmp", outImg);

    //start = chrono::high_resolution_clock::now();
    //SIFTGPUMatcher.MatchGuided(keypoints1, keypoints2, descriptors1, descriptors2, gpuTwoViewGeometry);
    //end = chrono::high_resolution_clock::now();
    //cout << "GPU guide match时间: " << chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << endl;
    //cout << "GPU指导匹配后的数量: " << gpuTwoViewGeometry.inlierMatches.size() << endl;
    //vector<cv::DMatch> dmatches6;
    //for (size_t i = 0; i < gpuTwoViewGeometry.inlierMatches.size(); i++)
    //{
    //    dmatches6.push_back(cv::DMatch(gpuTwoViewGeometry.inlierMatches[i].point2DIndex1, gpuTwoViewGeometry.inlierMatches[i].point2DIndex2, 0));
    //}
    //cv::drawMatches(img1, keypointsCV1, img2, keypointsCV2, dmatches6, outImg);
    //cv::imwrite("11GuideMatch_GPU.bmp", outImg);










	return 0;
}









