#include <iostream>
#include <string>
#include <vector>
#include <io.h>
#include <windows.h>
#include <chrono>
#include <thread>
#include <Windows.h>


#include <opencv2/opencv.hpp>
#include <thirdparty/SIFT_GPU/SiftGPU.h>
#include <Eigen/Core>
#include <GL/gl.h>
#include <src/Base/Base.h>
#include <opencv2/sfm.hpp>

using namespace std;
class CSiftExtractorOptions
{
public:
    int MaxThread = -1; // Number of threads for feature extraction.
    bool IsUseGPU = true;
    string GPUIndex = "-1";
    int MaxImageSize = 3200;
    int MaxImageFeature = 8192;
    int FirstOctave = -1;                            // First octave in the pyramid, i.e. -1 upsamples the image by one level.
    int OctavesNum = 4;
    int OctavesResolution = 3;                       // Number of levels per octave.
    double PeakThreshold = 0.02 / OctavesResolution; // 检测的峰值阈值
    double EdgeThreshold = 10.0;                     // 检测的边缘阈值
    bool IsEstimateAffine = false;                   // 估计SIFT特征的仿射形状, 采用有方向的椭圆形式, 而不是原始的SIFT方法中估计有方向的圆盘形状
    int MaxOrientationsNum = 2;                      // 如果不进行IsEstimateAffine，则每个关键点的最大方向数量
    bool IsUpright = false;                          // 将竖直特征的方向修正为0
    bool IsDarkAdapt = false;                        // 是否根据图像的暗度调整特征检测(仅适用于OpenGL版本的SiftGPU)
    
    // 域大小池化参数. 域大小池化计算在检测到的尺度周围多个尺度上的平均SIFT描述符, 在表现上优于其他SIFT变种和学习到的描述符
    bool IsDomainSizePooling = false;
    double DspMinScale = 1.0 / 6.0;
    double DspMaxScale = 3.0;
    int DspScaleNum = 10;

    bool IsForceUseVLFeat = false; // 是否强制使用协变的VLFeat实现. 否则, 只有在启用IsEstimateAffine或IsDomainSizePooling时才会使用协变实现, 因为普通的Sift实现更快速.

    enum class CNormalization
    {
        L1_ROOT, // 对每个描述符进行L1归一化，然后进行逐元素平方根运算. 通常优于标准的L2归一化. 《Three things everyone should know to improve object retrieval》
        L2       // 每个向量都进行L2归一化
    };
    CNormalization Normalization = CNormalization::L1_ROOT;
};

struct CFeatureKeypoint
{
    float x, y, a11, a12, a21, a22;
    CFeatureKeypoint()
    {
        CFeatureKeypoint(0, 0);
    }
    CFeatureKeypoint(float x, float y)
    {
        this->x = x;
        this->y = y;
        a11 = 1;
        a12 = 0;
        a21 = 0;
        a22 = 1;
    }
    CFeatureKeypoint(float x, float y, float scale, float orientation)
    {
        this->x = x;
        this->y = y;
        const float scale_cos_orientation = scale * cos(orientation);
        const float scale_sin_orientation = scale * sin(orientation);
        a11 = scale_cos_orientation;
        a12 = -scale_sin_orientation;
        a21 = scale_sin_orientation;
        a22 = scale_cos_orientation;
    }
    CFeatureKeypoint(float x, float y, float a11, float a12, float a21, float a22)
    {
        this->x = x;
        this->y = y;
        this->a11 = a11;
        this->a12 = a12;
        this->a21 = a21;
        this->a22 = a22;
    }
};
typedef vector<CFeatureKeypoint> CFeatureKeypoints;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CFeatureDescriptors;

template <typename T1, typename T2>
T2 TruncateCast(const T1 value) 
{
    return static_cast<T2>(std::min(static_cast<T1>(std::numeric_limits<T2>::max()), std::max(static_cast<T1>(std::numeric_limits<T2>::min()), value)));
}


class CSiftGPUExtractor
{
public:
    explicit CSiftGPUExtractor(const CSiftExtractorOptions& options) : SiftExtractorOptions(options)
    {

    }
    static unique_ptr<CSiftGPUExtractor> Create(const CSiftExtractorOptions& options)
    {
        vector<int> GPUIndex = { 0 };
        vector<string> Args;
        Args.push_back("./sift_gpu");
        Args.push_back("-cuda");
        Args.push_back(to_string(GPUIndex[0]));
        Args.push_back("-v");
        Args.push_back("0");

        const int CompensationFactor = 1 << -min(0, options.FirstOctave);
        Args.push_back("-maxd");
        Args.push_back(to_string(options.MaxImageSize * CompensationFactor));

        // Keep the highest level features.
        Args.push_back("-tc2");
        Args.push_back(to_string(options.MaxImageFeature));

        // First octave level.
        Args.push_back("-fo");
        Args.push_back(to_string(options.FirstOctave));

        // Number of octave levels.
        Args.push_back("-d");
        Args.push_back(to_string(options.OctavesResolution));

        // Peak threshold.
        Args.push_back("-t");
        Args.push_back(to_string(options.PeakThreshold));

        // Edge threshold.
        Args.push_back("-e");
        Args.push_back(to_string(options.EdgeThreshold));

        Args.push_back("-mo");
        Args.push_back(to_string(options.MaxOrientationsNum));

        vector<const char*> SiftGPUArgs;
        SiftGPUArgs.reserve(Args.size());
        for (const auto& arg : Args) 
        {
            SiftGPUArgs.push_back(arg.c_str());
        }

        unique_ptr<CSiftGPUExtractor> extractor = make_unique<CSiftGPUExtractor>(options);

        // Note that the SiftGPU object is not movable (for whatever reason). If we instead create the object here and move it to the constructor, the program segfaults inside SiftGPU.
        extractor->SiftGPU.ParseParam(SiftGPUArgs.size(), SiftGPUArgs.data());

        extractor->SiftGPU.gpu_index = GPUIndex[0];
        if (extractor->SiftGPU.VerifyContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED)
        {
            return nullptr;
        }

        return extractor;
    }
    bool Extract(const cv::Mat& Image, CFeatureKeypoints& Keypoints, CFeatureDescriptors& Descriptors)
    {
        const int compensation_factor = 1 << -min(0, SiftExtractorOptions.FirstOctave);

        int scan_width = Image.step;
        int bpp = Image.elemSize() * 8;
        vector<uint8_t> raw_bits(Image.total() * Image.elemSize(), 0);
        memcpy(raw_bits.data(), Image.data, raw_bits.size());
        const int code = SiftGPU.RunSIFT(scan_width, Image.rows, raw_bits.data(), GL_LUMINANCE, GL_UNSIGNED_BYTE);
        if (code != 1)
        {
            return false;
        }
        const size_t FeatureNum = static_cast<size_t>(SiftGPU.GetFeatureNum());
        vector<SiftKeypoint> KeypointsBuffer(FeatureNum);
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DescriptorsBuffer(FeatureNum, 128);
        SiftGPU.GetFeatureVector(KeypointsBuffer.data(), DescriptorsBuffer.data());
        Keypoints.resize(FeatureNum);

        for (size_t i = 0; i < FeatureNum; i++)
        {
            Keypoints[i].x = KeypointsBuffer[i].x;
            Keypoints[i].y = KeypointsBuffer[i].y;

            float scale_cos_orientation = KeypointsBuffer[i].s * cos(KeypointsBuffer[i].o);
            float scale_sin_orientation = KeypointsBuffer[i].s * sin(KeypointsBuffer[i].o);
            Keypoints[i].a11 = scale_cos_orientation;
            Keypoints[i].a12 = -scale_sin_orientation;
            Keypoints[i].a21 = scale_sin_orientation;
            Keypoints[i].a22 = scale_cos_orientation;
        }
        
        for (Eigen::MatrixXf::Index r = 0; r < DescriptorsBuffer.rows(); r++)
        {
            DescriptorsBuffer.row(r) *= 1 / DescriptorsBuffer.row(r).lpNorm<1>();
            DescriptorsBuffer.row(r) = DescriptorsBuffer.row(r).array().sqrt();
        }

        Descriptors = CFeatureDescriptors(FeatureNum, 128);
        for (Eigen::MatrixXf::Index r = 0; r < DescriptorsBuffer.rows(); r++)
        {
            for (Eigen::MatrixXf::Index c = 0; c < DescriptorsBuffer.cols(); c++)
            {
                const float scaled_value = std::round(512.0f * DescriptorsBuffer(r, c));
                Descriptors(r, c) = TruncateCast<float, uint8_t>(scaled_value);
            }
        }
        return true;
    }

private:
    CSiftExtractorOptions SiftExtractorOptions;
    SiftGPU SiftGPU;
};

void Extract(string Path)
{
    cv::Mat Image_Gray = cv::imread(Path, cv::IMREAD_GRAYSCALE);
    if (Image_Gray.empty())
    {
        cerr << "无法读取图像文件" << endl;
        return;
    }

    CSiftExtractorOptions options;
    unique_ptr<CSiftGPUExtractor> Extractor(CSiftGPUExtractor::Create(options));
    CFeatureKeypoints Keypoints;
    CFeatureDescriptors Descriptors;
    Extractor->Extract(Image_Gray, Keypoints, Descriptors);

    cv::Mat Image = cv::imread(Path, cv::IMREAD_ANYCOLOR);
    for (const auto& point : Keypoints)
    {
        cv::circle(Image, cv::Point(point.x, point.y), 5, cv::Scalar(0, 0, 255), -1); // Red circle
    }
    string FileName = Path.substr(Path.find_last_of('/') + 1);
    cv::imwrite(FileName, Image);
}


int main(int argc, char* argv[])
{

    vector<string> ImagePaths;
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00001.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00002.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00003.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00004.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00005.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00006.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00007.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00008.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00009.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00010.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00011.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00012.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00013.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00014.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00015.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00016.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00017.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00018.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00019.JPG");
    ImagePaths.push_back("D:/Dataset/Selfmade/磨山/四人石像/DSC00020.JPG");

    thread t1(&Extract, ImagePaths[0]);
    thread t2(&Extract, ImagePaths[1]);
    thread t3(&Extract, ImagePaths[2]);
    thread t4(&Extract, ImagePaths[3]);
    thread t5(&Extract, ImagePaths[4]);
    /*thread t6(&Extract, ImagePaths[5]);
    thread t7(&Extract, ImagePaths[6]);
    thread t8(&Extract, ImagePaths[7]);
    thread t9(&Extract, ImagePaths[8]);
    thread t10(&Extract, ImagePaths[9]);*/
    /*thread t11(&Extract, ImagePaths[10]);
    thread t12(&Extract, ImagePaths[11]);
    thread t13(&Extract, ImagePaths[12]);
    thread t14(&Extract, ImagePaths[13]);
    thread t15(&Extract, ImagePaths[14]);*/
    /*thread t16(&Extract, ImagePaths[15]);
    thread t17(&Extract, ImagePaths[16]);
    thread t18(&Extract, ImagePaths[17]);
    thread t19(&Extract, ImagePaths[18]);
    thread t20(&Extract, ImagePaths[19]);*/
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    /*t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();*/
    /*t11.join();
    t12.join();
    t13.join();
    t14.join();
    t15.join();*/
    /*t16.join();
    t17.join();
    t18.join();
    t19.join();
    t20.join();*/

    return 0;
}










