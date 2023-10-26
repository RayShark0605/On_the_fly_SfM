#include "FeatureExtractor.h"
#include <GL/GL.h>
using namespace std;

CSIFTCPUExtractor::CSIFTCPUExtractor(const CSIFTExtractionOptions& options) :sift(nullptr, &vl_sift_delete)
{
	options.CheckOptions();
	Check(!options.isEstimateAffineShape);
	Check(!options.isDomainSizePooling);
	this->options = options;
}
bool CSIFTCPUExtractor::Extract(const string& imagePath, CKeypoints& keypoints, CSIFTDescriptors& descriptors)
{
#ifdef SIFT_Use_OpenCV
    Check(IsFileExists(imagePath));
    cv::Mat image = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
    if (image.empty())
    {
        return false;
    }
    return Extract_OpenCV(image, keypoints, descriptors);

#else
    size_t originWidth, originHeight, currentWidth, currentHeight;
    vector<uchar> imageData = ReadImage(imagePath, options.maxImageSize, originWidth, originHeight, currentWidth, currentHeight);
    if (!Extract_VLFeat(imageData, currentWidth, currentHeight, keypoints, descriptors))
    {
        return false;
    }
    ScaleKeypoints(keypoints, currentWidth, currentHeight, originWidth, originHeight);
    return true;
#endif
}
bool CSIFTCPUExtractor::Extract_VLFeat(const std::vector<uchar>& imageData, size_t width, size_t height, CKeypoints& keypoints, CSIFTDescriptors& descriptors)
{
    Check(imageData.size() == width * height); // 确保是灰度影像数据

    // 如果sift对象不存在或尺寸不对, 就重新创建一个
    if (!sift || sift->width != width || sift->height != height)
    {
        sift = VlSiftType(vl_sift_new(width, height, options.numOctaves, options.octaveResolution, -1), &vl_sift_delete);
        if (!sift)
        {
            return false; // 创建失败, 返回false
        }
    }


    // 设置sift的峰值和边缘阈值
    vl_sift_set_peak_thresh(sift.get(), options.peakThreshold);
    vl_sift_set_edge_thresh(sift.get(), options.edgeThreshold);

    vector<size_t> numFeaturesPerLevel;          // 每个DOG层中检测到的关键点数量
    vector<CKeypoints> keypointsPerLevel;         // 每个DOG层中检测到的关键点信息
    vector<CSIFTDescriptors> descriptorsPerLevel; // 每个DOG层中关键点对应的SIFT描述子
    bool isFirstOctave = true; // 标记是否是第一个尺度(octave)

    // 循环处理每个尺度(octave)
    while (true)
    {
        if (isFirstOctave) // 针对第一个尺度(octave)进行特殊处理
        {
            // 将图像数据转换为浮点数, 并归一化到[0, 1]
            vector<float> imageData_Float(imageData.size());
            for (size_t i = 0; i < imageData.size(); i++)
            {
                imageData_Float[i] = imageData[i] * 1.0 / 255.0f;
            }

            // 处理第一个尺度(octave)
            if (vl_sift_process_first_octave(sift.get(), imageData_Float.data()))
            {
                break;
            }
            isFirstOctave = false;
        }
        else // 处理接下来的尺度(octave)
        {
            if (vl_sift_process_next_octave(sift.get()))
            {
                break; // 如果处理失败, 跳出循环
            }
        }

        // 在当前尺度(octave)下检测关键点
        vl_sift_detect(sift.get());

        // 获取检测到的关键点
        const VlSiftKeypoint* keypoints_VLFeat = vl_sift_get_keypoints(sift.get());
        const int numKeyPoints = vl_sift_get_nkeypoints(sift.get());
        if (numKeyPoints == 0) // 如果没有检测到关键点, 则继续下一个尺度(octave)
        {
            continue;
        }

        // 对不同的DOG(Difference of Gaussians)层进行特征提取
        size_t levelIndex = 0; // 当前DOG层的关键点计数
        int preLevel = -1;  // 上一个遍历的DOG层
        CSIFTDescriptors_Float desc(1, 128); // SIFT描述子

        // 遍历当前尺度下的所有关键点
        for (int i = 0; i < numKeyPoints; i++)
        {
            if (keypoints_VLFeat[i].is != preLevel) // 检查是否进入了新的DOG层
            {
                if (i > 0)
                {
                    // 调整前一个DOG层的容器大小
                    keypointsPerLevel.back().resize(levelIndex);
                    descriptorsPerLevel.back().conservativeResize(levelIndex, 128);
                }

                // 为新的DOG层添加容器
                levelIndex = 0;
                numFeaturesPerLevel.push_back(0);
                keypointsPerLevel.emplace_back(options.maxNumOrientations * numKeyPoints);
                descriptorsPerLevel.emplace_back(options.maxNumOrientations * numKeyPoints, 128);
            }

            // 更新当前DOG层的关键点数量
            numFeaturesPerLevel.back()++;
            preLevel = keypoints_VLFeat[i].is;  // 更新前一个DOG层的标记

            // 提取特征方向
            double angles[4];      // 存储方向的数组
            int numOrientations;  // 方向数量
            if (options.isUpRight) // 如果选项设置为"isUpRight", 只使用一个方向
            {
                numOrientations = 1;
                angles[0] = 0.0;
            }
            else  // 否则, 计算多个方向
            {
                numOrientations = vl_sift_calc_keypoint_orientations(sift.get(), angles, &keypoints_VLFeat[i]);
            }

            // 选择要使用的方向数量. 这里与SiftGPU不同, SiftGPU选择全局最大值作为方向, 而此处选择前两个局部最大值作为方向. 目前尚不清楚哪种方法更好
            const int numUsedOrientations = min((size_t)numOrientations, options.maxNumOrientations);

            // 针对每个选定的方向, 提取描述子
            for (int orientationIndex = 0; orientationIndex < numUsedOrientations; orientationIndex++)
            {
                // 创建并存储新的关键点对象
                keypointsPerLevel.back()[levelIndex] = CKeypoint(keypoints_VLFeat[i].x + 0.5f, keypoints_VLFeat[i].y + 0.5f, keypoints_VLFeat[i].sigma, angles[orientationIndex]);

                // 计算并存储新的描述子
                vl_sift_calc_keypoint_descriptor(sift.get(), desc.data(), &keypoints_VLFeat[i], angles[orientationIndex]);

                // 根据选项进行描述子的归一化处理
                if (options.normalizationType == CSIFTNormalizationType::L2)
                {
                    L2Normalize(desc);
                }
                else if (options.normalizationType == CSIFTNormalizationType::L1_ROOT)
                {
                    L1RootNormalize(desc);
                }
                else
                {
                    Check(false, "Normalization type not supported");
                }

                // 转换描述子到无符号字节格式并存储
                descriptorsPerLevel.back().row(levelIndex) = SIFTDescriptorsFloatToUnsignedChar(desc);

                // 更新当前DOG层的关键点计数
                levelIndex++;
            }
        }

        // 调整最后一个DOG层的容器大小
        keypointsPerLevel.back().resize(levelIndex);
        descriptorsPerLevel.back().conservativeResize(levelIndex, 128);
    }

    int firstLevelIndexToKeep = 0;       // 第一个需要保留的DOG层的索引
    int numFeatures = 0;                 // 关键点的总数量
    int numFeaturesWithOrientations = 0; // 包括多个方向的关键点的总数量

    // 倒序遍历所有DOG层, 以确定哪些关键点将被保留
    for (int i = keypointsPerLevel.size() - 1; i >= 0; i--)
    {
        numFeatures += numFeaturesPerLevel[i]; // 更新关键点总数
        numFeaturesWithOrientations += keypointsPerLevel[i].size(); // 更新包括多个方向的关键点总数

        // 如果关键点总数超过了设定的最大数量，记录下第一个需要保留的DOG层
        if (numFeatures > options.maxNumFeatures)
        {
            firstLevelIndexToKeep = i;
            break;
        }
    }

    // 提取要保留的关键点
    size_t currentIndex = 0;
    keypoints.resize(numFeaturesWithOrientations); // 调整输出关键点数组的大小

    // 遍历从第一个需要保留的DOG层开始的所有DOG层
    for (size_t i = firstLevelIndexToKeep; i < keypointsPerLevel.size(); i++)
    {
        for (size_t j = 0; j < keypointsPerLevel[i].size(); j++)
        {
            keypoints[currentIndex] = keypointsPerLevel[i][j]; // 保存关键点
            currentIndex++; // 更新输出数组的索引
        }
    }

    // 计算检测到的关键点的描述子
    size_t k = 0;
    descriptors.resize(numFeaturesWithOrientations, 128); // 调整输出描述子矩阵的大小

    // 同样从第一个需要保留的DOG层开始，遍历所有DOG层
    for (size_t i = firstLevelIndexToKeep; i < keypointsPerLevel.size(); i++)
    {
        for (size_t j = 0; j < keypointsPerLevel[i].size(); j++)
        {
            descriptors.row(k) = descriptorsPerLevel[i].row(j); // 保存描述子
            k++; // 更新输出矩阵的索引
        }
    }
    descriptors = VLFeatToOriginFormat(descriptors); // 转换VLFeat格式的描述子到UBC格式

    ExtractTopScaleFeatures(keypoints, descriptors, options.maxNumFeatures);
    return true;
}
bool CSIFTCPUExtractor::Extract_OpenCV(const cv::Mat& image, CKeypoints& keypoints_Output, CSIFTDescriptors& descriptors_Output)
{
    cv::Mat descriptors;
    vector<cv::KeyPoint> keypoints;
    cv::Ptr<cv::SIFT> sift = cv::SIFT::create(options.maxNumFeatures, options.numOctaves, options.peakThreshold, options.edgeThreshold);
    sift->detectAndCompute(image, cv::noArray(), keypoints, descriptors);
    if (keypoints.empty())
    {
        return false;
    }
    size_t numFeatures = keypoints.size();
    Check(descriptors.type() == CV_32F && descriptors.cols == 128 && descriptors.rows == numFeatures, "Invalid input descriptors. Must be of type CV_32F and have 128 columns.");
    
    // 将OpenCV格式的特征点和描述子格式转换成本程序的格式
    keypoints_Output.resize(numFeatures);
    descriptors_Output.resize(numFeatures, 128);
    for (size_t i = 0; i < numFeatures; i++)
    {
        keypoints_Output[i] = move(CKeypoint(keypoints[i].pt.x, keypoints[i].pt.y, keypoints[i].size, keypoints[i].angle));

        cv::Mat row = descriptors.row(i);
        if (options.normalizationType == CSIFTNormalizationType::L1_ROOT) // 做L1_ROOT归一化
        {
            const double l1norm = cv::norm(row, cv::NORM_L1);
            row /= l1norm;
            for (int c = 0; c < 128; c++)
            {
                row.at<float>(0, c) = sqrt(row.at<float>(0, c));
            }
        }
        else if (options.normalizationType == CSIFTNormalizationType::L2) // 做L2归一化
        {
            const double l2norm = cv::norm(row, cv::NORM_L2);
            row /= l2norm;
        }
        else
        {
            Check(false, "Normalization type not supported");
        }

        for (int c = 0; c < 128; c++)
        {
            float originValue = row.at<float>(0, c);
            Check(originValue >= 0 && originValue <= 0.5);
            const float scaledValue = round(512.0f * originValue);
            descriptors_Output(i, c) = TruncateCast<float, uint8_t>(scaledValue);
        }
    }
    return true;
}

CSIFTGPUExtractor::CSIFTGPUExtractor(const CSIFTExtractionOptions& options)
{
    options.CheckOptions();
    Check(!options.isEstimateAffineShape);
    Check(!options.isDomainSizePooling);
    this->options = options;

    vector<const char*> siftGPU_Args;
    siftGPU_Args.push_back("./sift_gpu");
    siftGPU_Args.push_back("-cuda");

    size_t cudaDeviceIndex = SetBestCudaDevice();

    string cudaDeviceIndex_String = to_string(cudaDeviceIndex);
    string maxImageSize_String = to_string(options.maxImageSize * 2);
    string maxNumFeatures_String = to_string(options.maxNumFeatures);
    string octaveResolution_String = to_string(options.octaveResolution);
    string peakThreshold_String = to_string(options.peakThreshold);
    string edgeThreshold_String = to_string(options.edgeThreshold);
    string maxNumOrientations_String = to_string(options.maxNumOrientations);


    siftGPU_Args.push_back(cudaDeviceIndex_String.c_str());

    siftGPU_Args.push_back("-v");
    siftGPU_Args.push_back("0");
    siftGPU_Args.push_back("-maxd");
    siftGPU_Args.push_back(maxImageSize_String.c_str());
    siftGPU_Args.push_back("-tc2");
    siftGPU_Args.push_back(maxNumFeatures_String.c_str());
    siftGPU_Args.push_back("-fo");
    siftGPU_Args.push_back("-1");
    siftGPU_Args.push_back("-d");
    siftGPU_Args.push_back(octaveResolution_String.c_str());
    siftGPU_Args.push_back("-t");
    siftGPU_Args.push_back(peakThreshold_String.c_str());
    siftGPU_Args.push_back("-e");
    siftGPU_Args.push_back(edgeThreshold_String.c_str());

    if (options.isUpRight)
    {
        siftGPU_Args.push_back("-ofix");
        siftGPU_Args.push_back("-mo");
        siftGPU_Args.push_back("1");
    }
    else
    {
        siftGPU_Args.push_back("-mo");
        siftGPU_Args.push_back(maxNumOrientations_String.c_str());
    }
    siftGPU.ParseParam(siftGPU_Args.size(), siftGPU_Args.data());
    siftGPU.gpu_index = cudaDeviceIndex;

    Check(siftGPU.VerifyContextGL() == SiftGPU::SIFTGPU_FULL_SUPPORTED);
}
bool CSIFTGPUExtractor::Extract(const std::string& imagePath, CKeypoints& keypoints, CSIFTDescriptors& descriptors)
{
    size_t originWidth, originHeight, currentWidth, currentHeight;
    vector<uchar> imageData = ReadImage(imagePath, options.maxImageSize, originWidth, originHeight, currentWidth, currentHeight);
    if (!Extract_SIFTGPU(imageData, currentWidth, currentHeight, keypoints, descriptors))
    {
        return false;
    }
    ScaleKeypoints(keypoints, currentWidth, currentHeight, originWidth, originHeight);
    return true;
}
bool CSIFTGPUExtractor::Extract_SIFTGPU(const std::vector<uchar>& imageData, size_t width, size_t height, CKeypoints& keypoints, CSIFTDescriptors& descriptors)
{
    Check(imageData.size() == width * height); // 确保是灰度影像数据
    Check(options.maxImageSize * 2 == siftGPU.GetMaxDimension());

    const int code = siftGPU.RunSIFT(width, height, imageData.data(), GL_LUMINANCE, GL_UNSIGNED_BYTE);
    if (code != 1)
    {
        return false;
    }

    const size_t numFeatures = siftGPU.GetFeatureNum();
    vector<SiftKeypoint> keypointsBuffer(numFeatures);
    CSIFTDescriptors_Float descriptorsFloat(numFeatures, 128);
    siftGPU.GetFeatureVector(keypointsBuffer.data(), descriptorsFloat.data());

    keypoints.resize(numFeatures);
    for (size_t i = 0; i < numFeatures; i++)
    {
        keypoints[i] = move(CKeypoint(keypointsBuffer[i].x, keypointsBuffer[i].y, keypointsBuffer[i].s, keypointsBuffer[i].o));
    }
    if (options.normalizationType == CSIFTNormalizationType::L1_ROOT)
    {
        L1RootNormalize(descriptorsFloat);
    }
    else if (options.normalizationType == CSIFTNormalizationType::L2)
    {
        L2Normalize(descriptorsFloat);
    }
    else
    {
        Check(false, "Normalization type not supported");
    }
    descriptors = SIFTDescriptorsFloatToUnsignedChar(descriptorsFloat);
    ExtractTopScaleFeatures(keypoints, descriptors, options.maxNumFeatures);
    return true;
}