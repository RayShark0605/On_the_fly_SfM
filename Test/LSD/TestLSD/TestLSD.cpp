#include <iostream>
#include <string>
#include <vector>
#include <io.h>
#include <windows.h>
#include <chrono>


#include <opencv2/opencv.hpp>
#include <thirdparty/LSD/lsd.h>

using namespace std;

vector<string> findImageFiles(const string& folderPath) 
{
    vector<string> imagePaths;
    struct _finddata_t fileData;
    intptr_t handle = _findfirst((folderPath + "/*.*").c_str(), &fileData);

    if (handle != -1) 
    {
        do 
        {
            string fileName = fileData.name;
            if (fileName != "." && fileName != "..") 
            {
                string filePath = folderPath + "/" + fileName;

                // 使用OpenCV来判断文件是否为影像文件
                cv::Mat image = cv::imread(filePath);
                if (!image.empty()) 
                {
                    imagePaths.push_back(filePath);
                }
            }
        } while (_findnext(handle, &fileData) == 0);

        _findclose(handle);
    }
    return imagePaths;
}



int main(int argc, char* argv[])
{
    string TestImagePath = "D:/Dataset/TUM/Structure_vs._Texture/rgbd_dataset_freiburg3_nostructure_texture_far/rgb/1341839847.173932.png";
    cv::Mat Image_Gray = cv::imread(TestImagePath, cv::IMREAD_GRAYSCALE);
    if (Image_Gray.empty())
    {
        cerr << "无法读取图像文件" << endl;
        return -1;
    }
    vector<uint8_t> imageData(Image_Gray.begin<uint8_t>(), Image_Gray.end<uint8_t>());
    vector<double> bitmap_data_double(imageData.begin(), imageData.end());
    int LineNum = 0;
    double* Result = lsd(&LineNum, bitmap_data_double.data(), Image_Gray.cols, Image_Gray.rows);

    cv::Mat Image = cv::imread(TestImagePath, cv::IMREAD_ANYCOLOR);
    for (int i = 0; i < LineNum; i++)
    {
        double x1 = Result[7 * i];
        double y1 = Result[7 * i + 1];
        double x2 = Result[7 * i + 2];
        double y2 = Result[7 * i + 3];
        double width = Result[7 * i + 4];

        cv::line(Image, cv::Point(cvRound(x1), cvRound(y1)), cv::Point(cvRound(x2), cvRound(y2)), cv::Scalar(0, 0, 255), 1);
    }

    cv::Ptr<cv::LineSegmentDetector> lsd_Detector = cv::createLineSegmentDetector(cv::LSD_REFINE_STD);
    vector<cv::Vec4f> lines;
    lsd_Detector->detect(Image_Gray, lines);
    cv::Mat Image_Color = cv::imread(TestImagePath, cv::IMREAD_ANYCOLOR);
    cv::Mat drawnLines(Image_Color);
    lsd_Detector->drawSegments(drawnLines, lines);




    cv::imshow("LSD Line Detection_LSD", Image);
    cv::imshow("LSD Line Detection_OpenCV", Image_Color);
    cv::waitKey(0);
    free(Result);


    cout << "测试与OpenCV的速度对比!" << endl;
    string ImageDir = "D:/Dataset/TUM/Structure_vs._Texture/rgbd_dataset_freiburg3_nostructure_texture_far/rgb";
    vector<string> imagePaths = findImageFiles(ImageDir);

    cout << "LSD时间: " << endl;
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < imagePaths.size(); i++)
    {
        Image_Gray = cv::imread(imagePaths[i], cv::IMREAD_GRAYSCALE);
        if (Image_Gray.empty())
        {
            cerr << "无法读取图像文件" << endl;
            return -1;
        }
        vector<uint8_t> imageData(Image_Gray.begin<uint8_t>(), Image_Gray.end<uint8_t>());
        vector<double> bitmap_data_double(imageData.begin(), imageData.end());
        int LineNum = 0;
        Result = lsd(&LineNum, bitmap_data_double.data(), Image_Gray.cols, Image_Gray.rows);
        cout << LineNum << endl;
        free(Result);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "程序运行时间：" << duration.count() << " 毫秒" << endl;

    cout << "OpenCV时间: " << endl;
    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < imagePaths.size(); i++)
    {
        cv::Mat Image_Gray = cv::imread(imagePaths[i], cv::IMREAD_GRAYSCALE);
        cv::Ptr<cv::LineSegmentDetector> lsd = cv::createLineSegmentDetector(cv::LSD_REFINE_STD);
        vector<cv::Vec4f> lines;
        lsd->detect(Image_Gray, lines);
        cout << lines.size() << endl;
    }
    end = std::chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "程序运行时间：" << duration.count() << " 毫秒" << endl;

    return 0;
}


