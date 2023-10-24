#include <vector>
#include <string>
#include <iostream>
#include <tbb/concurrent_queue.h>
#include "../../src/Base/Log.h"
#include "../../src/Workflow/exif.h"
#include "../../src/Workflow/ThreadPool.h"
#include "../../src/Workflow/Workflow.h"
#include <QtWidgets/QApplication>
#include "../../src/UI/ModelViewerWidget.h"

using namespace std;

COptions options;
CRAMDatabase database;
void Func()
{
	this_thread::sleep_for(chrono::seconds(5));
	cout << "Æô¶¯!" << endl;
	string imageDirPath = "C:/Users/XiaRui/Desktop/PTRS_Project/ColmapDebugImages2";
	const vector<string> images = GetImageFileList(imageDirPath);

	const size_t num = 2;
	for (size_t i = 0; i < num; i++)
	{
		const string imagePath = images[i];
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
			CWorkflow::EstimateTwoViewGeometry(database, i, imageID, options.twoViewGeometryOptions);
		}
		auto end = chrono::high_resolution_clock::now();
		size_t timeConsuming = chrono::duration_cast<chrono::milliseconds>(end - start).count();
		const size_t numMatchedImages = database.GetTwoViewGeometryImagesNum(imageID);
		CLog::Log((boost::format("[%1%ms] %2% have %3% matched images!") % timeConsuming % imageName % numMatchedImages).str());
	}
}


int main(int argc, char* argv[])
{
	CLog log;
	QApplication app(argc, argv);

	options.twoViewGeometryOptions.isComputeRelativePose = true;
	CModelViewerWidget modelViewerWidget(nullptr, &options, &database);
	modelViewerWidget.show();

	thread t(Func);
	t.join();

	//const CTwoViewGeometry& twoViewGeometry = database.GetTwoViewGeometry(0, 1);
	//database.GetImage(1).GetWorldToCamera() = twoViewGeometry.image1ToImage2;
	database.GetImage(1).GetWorldToCamera(0) = CRigid3D(Eigen::Quaterniond::Identity(), Eigen::Vector3d(10, 10, 10));
	modelViewerWidget.ShowImages({ 0,1 }, 0);


	return app.exec();
}








