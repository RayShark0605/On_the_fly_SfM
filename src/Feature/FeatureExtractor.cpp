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
    Check(imageData.size() == width * height); // ȷ���ǻҶ�Ӱ������

    // ���sift���󲻴��ڻ�ߴ粻��, �����´���һ��
    if (!sift || sift->width != width || sift->height != height)
    {
        sift = VlSiftType(vl_sift_new(width, height, options.numOctaves, options.octaveResolution, -1), &vl_sift_delete);
        if (!sift)
        {
            return false; // ����ʧ��, ����false
        }
    }


    // ����sift�ķ�ֵ�ͱ�Ե��ֵ
    vl_sift_set_peak_thresh(sift.get(), options.peakThreshold);
    vl_sift_set_edge_thresh(sift.get(), options.edgeThreshold);

    vector<size_t> numFeaturesPerLevel;          // ÿ��DOG���м�⵽�Ĺؼ�������
    vector<CKeypoints> keypointsPerLevel;         // ÿ��DOG���м�⵽�Ĺؼ�����Ϣ
    vector<CSIFTDescriptors> descriptorsPerLevel; // ÿ��DOG���йؼ����Ӧ��SIFT������
    bool isFirstOctave = true; // ����Ƿ��ǵ�һ���߶�(octave)

    // ѭ������ÿ���߶�(octave)
    while (true)
    {
        if (isFirstOctave) // ��Ե�һ���߶�(octave)�������⴦��
        {
            // ��ͼ������ת��Ϊ������, ����һ����[0, 1]
            vector<float> imageData_Float(imageData.size());
            for (size_t i = 0; i < imageData.size(); i++)
            {
                imageData_Float[i] = imageData[i] * 1.0 / 255.0f;
            }

            // �����һ���߶�(octave)
            if (vl_sift_process_first_octave(sift.get(), imageData_Float.data()))
            {
                break;
            }
            isFirstOctave = false;
        }
        else // ����������ĳ߶�(octave)
        {
            if (vl_sift_process_next_octave(sift.get()))
            {
                break; // �������ʧ��, ����ѭ��
            }
        }

        // �ڵ�ǰ�߶�(octave)�¼��ؼ���
        vl_sift_detect(sift.get());

        // ��ȡ��⵽�Ĺؼ���
        const VlSiftKeypoint* keypoints_VLFeat = vl_sift_get_keypoints(sift.get());
        const int numKeyPoints = vl_sift_get_nkeypoints(sift.get());
        if (numKeyPoints == 0) // ���û�м�⵽�ؼ���, �������һ���߶�(octave)
        {
            continue;
        }

        // �Բ�ͬ��DOG(Difference of Gaussians)�����������ȡ
        size_t levelIndex = 0; // ��ǰDOG��Ĺؼ������
        int preLevel = -1;  // ��һ��������DOG��
        CSIFTDescriptors_Float desc(1, 128); // SIFT������

        // ������ǰ�߶��µ����йؼ���
        for (int i = 0; i < numKeyPoints; i++)
        {
            if (keypoints_VLFeat[i].is != preLevel) // ����Ƿ�������µ�DOG��
            {
                if (i > 0)
                {
                    // ����ǰһ��DOG���������С
                    keypointsPerLevel.back().resize(levelIndex);
                    descriptorsPerLevel.back().conservativeResize(levelIndex, 128);
                }

                // Ϊ�µ�DOG���������
                levelIndex = 0;
                numFeaturesPerLevel.push_back(0);
                keypointsPerLevel.emplace_back(options.maxNumOrientations * numKeyPoints);
                descriptorsPerLevel.emplace_back(options.maxNumOrientations * numKeyPoints, 128);
            }

            // ���µ�ǰDOG��Ĺؼ�������
            numFeaturesPerLevel.back()++;
            preLevel = keypoints_VLFeat[i].is;  // ����ǰһ��DOG��ı��

            // ��ȡ��������
            double angles[4];      // �洢���������
            int numOrientations;  // ��������
            if (options.isUpRight) // ���ѡ������Ϊ"isUpRight", ֻʹ��һ������
            {
                numOrientations = 1;
                angles[0] = 0.0;
            }
            else  // ����, ����������
            {
                numOrientations = vl_sift_calc_keypoint_orientations(sift.get(), angles, &keypoints_VLFeat[i]);
            }

            // ѡ��Ҫʹ�õķ�������. ������SiftGPU��ͬ, SiftGPUѡ��ȫ�����ֵ��Ϊ����, ���˴�ѡ��ǰ�����ֲ����ֵ��Ϊ����. Ŀǰ�в�������ַ�������
            const int numUsedOrientations = min((size_t)numOrientations, options.maxNumOrientations);

            // ���ÿ��ѡ���ķ���, ��ȡ������
            for (int orientationIndex = 0; orientationIndex < numUsedOrientations; orientationIndex++)
            {
                // �������洢�µĹؼ������
                keypointsPerLevel.back()[levelIndex] = CKeypoint(keypoints_VLFeat[i].x + 0.5f, keypoints_VLFeat[i].y + 0.5f, keypoints_VLFeat[i].sigma, angles[orientationIndex]);

                // ���㲢�洢�µ�������
                vl_sift_calc_keypoint_descriptor(sift.get(), desc.data(), &keypoints_VLFeat[i], angles[orientationIndex]);

                // ����ѡ����������ӵĹ�һ������
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

                // ת�������ӵ��޷����ֽڸ�ʽ���洢
                descriptorsPerLevel.back().row(levelIndex) = SIFTDescriptorsFloatToUnsignedChar(desc);

                // ���µ�ǰDOG��Ĺؼ������
                levelIndex++;
            }
        }

        // �������һ��DOG���������С
        keypointsPerLevel.back().resize(levelIndex);
        descriptorsPerLevel.back().conservativeResize(levelIndex, 128);
    }

    int firstLevelIndexToKeep = 0;       // ��һ����Ҫ������DOG�������
    int numFeatures = 0;                 // �ؼ����������
    int numFeaturesWithOrientations = 0; // �����������Ĺؼ����������

    // �����������DOG��, ��ȷ����Щ�ؼ��㽫������
    for (int i = keypointsPerLevel.size() - 1; i >= 0; i--)
    {
        numFeatures += numFeaturesPerLevel[i]; // ���¹ؼ�������
        numFeaturesWithOrientations += keypointsPerLevel[i].size(); // ���°����������Ĺؼ�������

        // ����ؼ��������������趨�������������¼�µ�һ����Ҫ������DOG��
        if (numFeatures > options.maxNumFeatures)
        {
            firstLevelIndexToKeep = i;
            break;
        }
    }

    // ��ȡҪ�����Ĺؼ���
    size_t currentIndex = 0;
    keypoints.resize(numFeaturesWithOrientations); // ��������ؼ�������Ĵ�С

    // �����ӵ�һ����Ҫ������DOG�㿪ʼ������DOG��
    for (size_t i = firstLevelIndexToKeep; i < keypointsPerLevel.size(); i++)
    {
        for (size_t j = 0; j < keypointsPerLevel[i].size(); j++)
        {
            keypoints[currentIndex] = keypointsPerLevel[i][j]; // ����ؼ���
            currentIndex++; // ����������������
        }
    }

    // �����⵽�Ĺؼ����������
    size_t k = 0;
    descriptors.resize(numFeaturesWithOrientations, 128); // ������������Ӿ���Ĵ�С

    // ͬ���ӵ�һ����Ҫ������DOG�㿪ʼ����������DOG��
    for (size_t i = firstLevelIndexToKeep; i < keypointsPerLevel.size(); i++)
    {
        for (size_t j = 0; j < keypointsPerLevel[i].size(); j++)
        {
            descriptors.row(k) = descriptorsPerLevel[i].row(j); // ����������
            k++; // ����������������
        }
    }
    descriptors = VLFeatToOriginFormat(descriptors); // ת��VLFeat��ʽ�������ӵ�UBC��ʽ

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
    
    // ��OpenCV��ʽ��������������Ӹ�ʽת���ɱ�����ĸ�ʽ
    keypoints_Output.resize(numFeatures);
    descriptors_Output.resize(numFeatures, 128);
    for (size_t i = 0; i < numFeatures; i++)
    {
        keypoints_Output[i] = move(CKeypoint(keypoints[i].pt.x, keypoints[i].pt.y, keypoints[i].size, keypoints[i].angle));

        cv::Mat row = descriptors.row(i);
        if (options.normalizationType == CSIFTNormalizationType::L1_ROOT) // ��L1_ROOT��һ��
        {
            const double l1norm = cv::norm(row, cv::NORM_L1);
            row /= l1norm;
            for (int c = 0; c < 128; c++)
            {
                row.at<float>(0, c) = sqrt(row.at<float>(0, c));
            }
        }
        else if (options.normalizationType == CSIFTNormalizationType::L2) // ��L2��һ��
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
    Check(imageData.size() == width * height); // ȷ���ǻҶ�Ӱ������
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