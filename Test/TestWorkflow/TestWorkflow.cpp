#include <vector>
#include <string>
#include <iostream>
#include <tbb/concurrent_queue.h>
#include "../../src/Base/Log.h"
#include "../../src/Workflow/exif.h"
#include "../../src/Workflow/ThreadPool.h"
#include "../../src/Workflow/Workflow.h"

using namespace std;

vector<string> images;
tbb::concurrent_queue<string> q;
atomic_bool isStop;
CRAMDatabase database;
COptions options;
tbb::concurrent_queue<pair<size_t, size_t>> two;
atomic_size_t atomicCount;

void ProducerFunc()
{
	for (size_t i = 0; i < images.size(); i++)
	{
		this_thread::sleep_for(chrono::milliseconds(500));
		q.push(images[i]);
	}
}
void ConsumerFunc2()
{
	constexpr size_t invalidCount = numeric_limits<size_t>::max();
	while (!isStop)
	{
		this_thread::sleep_for(chrono::milliseconds(10));
		if (two.empty())
		{
			continue;
		}
		pair<size_t, size_t> matchPair({ invalidCount ,invalidCount });
		if (!two.try_pop(matchPair) || matchPair.first == invalidCount || matchPair.second == invalidCount)
		{
			continue;
		}
		CWorkflow::EstimateTwoViewGeometry(database, matchPair.first, matchPair.second, options.twoViewGeometryOptions);

		atomicCount--;
	}
}
void ConsumerFunc()
{
	while (!isStop)
	{
		this_thread::sleep_for(chrono::milliseconds(50));
		if (images.empty())
		{
			continue;
		}

		string imagePath = "";
		if (!q.try_pop(imagePath) || imagePath.empty())
		{
			continue;
		}
		const string imageName = GetFileName(imagePath);
		CLog::Log((boost::format("\nprocessing %1% ...") % imageName).str());
		
		auto start = chrono::high_resolution_clock::now();
		CWorkflow::ReadImage(imagePath, database);
		CWorkflow::SIFTExtract(imagePath, options.SIFTExtractionOptions, database);
		
		const size_t imageID = database.GetImageID(imageName);
		for (size_t i = 0; i < imageID; i++)
		{
			CWorkflow::SIFTMatch(database, i, imageID, options.SIFTMatchingOptions);
		}
		for (size_t i = 0; i < imageID; i++)
		{
			two.push({ i,imageID });
			atomicCount++;
		}
		while (!isStop && atomicCount > 0)
		{
			this_thread::sleep_for(chrono::milliseconds(100));
		}
		auto end = chrono::high_resolution_clock::now();
		size_t timeConsuming = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		const size_t numMatchedImages = database.GetTwoViewGeometryImagesNum(imageID);
		CLog::Log((boost::format("[%1%ms] %2% have %3% matched images!") % timeConsuming  % imageName % numMatchedImages).str());
	}
}


int main(int argc, char* argv[])
{
	CLog log;
	isStop = false;
	atomicCount = 0;

	cout << "Æô¶¯!" << endl;
	string imageDirPath = "C:/Users/XiaRui/Desktop/PTRS_Project/ColmapDebugImages2";
	images = GetImageFileList(imageDirPath);

	size_t numThread2 = 16;
	vector<thread> consumer2;
	consumer2.reserve(numThread2);
	for (size_t i = 0; i < numThread2; i++)
	{
		consumer2.emplace_back(ConsumerFunc2);
		consumer2.back().detach();
	}

	thread consumer(ConsumerFunc);
	consumer.detach();

	thread producer(ProducerFunc);
	producer.detach();


	this_thread::sleep_for(chrono::seconds(45));
	isStop = true;
	const vector<CImage>& databaseImages = database.GetAllImages();
	for (int i = 0; i < images.size(); i++)
	{
		const string leftImagePath = images[i];
		const string leftImageName = GetFileNameNoExtension(leftImagePath);
		const size_t leftImageID = database.GetImageID(GetFileName(leftImagePath));
		for (int j = i + 1; j < images.size(); j++)
		{
			const string rightImagePath = images[j];
			const string rightImageName = GetFileNameNoExtension(rightImagePath);
			const size_t rightImageID = database.GetImageID(GetFileName(rightImagePath));
			
			const size_t matchesNum = database.GetInlierMatchesNum(leftImageID, rightImageID);
			cout << leftImageName << "-" << rightImageName << ":" << matchesNum << endl;
			if (matchesNum > 0)
			{
				ExportMatches(leftImagePath, rightImagePath, database.GetImageKeypoints(leftImageID), database.GetImageKeypoints(rightImageID), database.GetTwoViewGeometry(leftImageID, rightImageID).inlierMatches, leftImageName + rightImageName + ".bmp");
			}
		}
	}

	/*string leftImagePath = "C:/Users/XiaRui/Desktop/PTRS_Project/images/DSC00001.JPG";
	string rightImagePath = "C:/Users/XiaRui/Desktop/PTRS_Project/images/DSC00026.JPG";
	size_t leftImageID = database.GetImageID(GetFileName(leftImagePath));
	size_t rightImageID = database.GetImageID(GetFileName(rightImagePath));
	cout << leftImageID << endl << rightImageID << endl;
	ExportMatches(leftImagePath, rightImagePath, database.GetImageKeypoints(leftImageID), database.GetImageKeypoints(rightImageID), database.GetTwoViewGeometry(leftImageID, rightImageID).inlierMatches, "0126.bmp");*/
	////ExportMatches(leftImagePath, rightImagePath, database.GetImageKeypoints(leftImageID), database.GetImageKeypoints(rightImageID), database.GetTwoViewGeometry(leftImageID, rightImageID).inlierMatches, "1516.bmp");
	return 0;

}












