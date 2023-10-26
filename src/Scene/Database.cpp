#include "Database.h"
#include <fstream>

using namespace std;

CRAMDatabase::CRAMDatabase() noexcept
{
	Clear();
	cameras.reserve(100);
	images.reserve(2000);
	matches.reserve(2000);
	twoViewGeometries.reserve(2000);
}
CRAMDatabase::~CRAMDatabase()
{
	Clear();
}
void CRAMDatabase::Clear() noexcept
{
	cameras.clear();
	images.clear();
	matches.clear();
	twoViewGeometries.clear();
}

size_t CRAMDatabase::GetCamerasNum() const noexcept
{
	return cameras.size();
}
bool CRAMDatabase::IsCameraExists(size_t cameraID) const
{
	return (cameras.find(cameraID) != cameras.end());
}
bool CRAMDatabase::IsCameraExists(const CCamera& camera) const
{
	for (const auto& pair : cameras)
	{
		if (pair.second == camera)
		{
			return true;
		}
	}
	return false;
}
size_t CRAMDatabase::GetCameraID(const CCamera& camera) const
{
	for (const auto& pair : cameras)
	{
		if (pair.second == camera)
		{
			return pair.first;
		}
	}
	Check(false, "This camera does not exist!");
}
const CCamera& CRAMDatabase::GetCamera(size_t cameraID) const
{
	auto it = cameras.find(cameraID);
	if (it != cameras.end())
	{
		return it->second;
	}
	Check(false, (boost::format("This camera does not exist! cameraID = %1%") % cameraID).str());
}
CCamera& CRAMDatabase::GetCamera(size_t cameraID)
{
	auto it = cameras.find(cameraID);
	if (it != cameras.end())
	{
		return it->second;
	}
	Check(false, (boost::format("This camera does not exist! cameraID = %1%") % cameraID).str());
}
vector<CCamera> CRAMDatabase::GetAllCameras() const
{
	vector<CCamera> result;
	result.reserve(cameras.size());
	for (const auto& pair : cameras)
	{
		result.push_back(pair.second);
	}
	return result;
}
size_t CRAMDatabase::AddCamera(const CCamera& camera)
{
	const size_t cameraID = cameras.size();
	cameras[cameraID] = camera;
	cameras[cameraID].SetCameraID(cameraID);
	return cameraID;
}

size_t CRAMDatabase::GetImagesNum() const noexcept
{
	return images.size();
}
bool CRAMDatabase::IsImageExists(size_t imageID) const
{
	return (images.find(imageID) != images.end());
}
bool CRAMDatabase::IsImageExists(const CImage& image) const
{
	for (const auto& pair : images)
	{
		if (pair.second == image)
		{
			return true;
		}
	}
	return false;
}
bool CRAMDatabase::IsImageExists(const string& imageName) const
{
	for (const auto& pair : images)
	{
		if (pair.second.GetImageName() == imageName)
		{
			return true;
		}
	}
	return false;
}
size_t CRAMDatabase::GetImageID(const string& imageName) const
{
	for (const auto& pair : images)
	{
		if (pair.second.GetImageName() == imageName)
		{
			return pair.first;
		}
	}
	return numeric_limits<size_t>::max();
}
const CCamera& CRAMDatabase::GetImageCamera(size_t imageID) const
{
	auto findImage = images.find(imageID);
	if (findImage != images.end())
	{
		const size_t cameraID = findImage->second.GetCameraID();
		auto findCamera = cameras.find(cameraID);
		if (findCamera != cameras.end())
		{
			return findCamera->second;
		}
		Check(false, (boost::format("This camera does not exist! cameraID = %1%") % cameraID).str());
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
const CImage& CRAMDatabase::GetImage(size_t imageID) const
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second;
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
CImage& CRAMDatabase::GetImage(size_t imageID)
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second;
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
vector<CImage> CRAMDatabase::GetAllImages() const
{
	vector<CImage> result;
	result.reserve(images.size());
	for (const auto& pair : images)
	{
		result.push_back(pair.second);
	}
	return result;
}
vector<string> CRAMDatabase::GetAllImageNames() const
{
	vector<string> result;
	result.reserve(images.size());
	for (const auto& pair : images)
	{
		result.push_back(pair.second.GetImageName());
	}
	return result;
}
size_t CRAMDatabase::AddImage(const CImage& image)
{
	const size_t imageID = images.size();
	images[imageID] = image;
	images[imageID].SetImageID(imageID);

	const size_t cameraID = images[imageID].GetCameraID();
	Check(cameras.find(cameraID) != cameras.end());
	images[imageID].Setup(cameras[cameraID]);

	return imageID;
}

size_t CRAMDatabase::GetImageKeypointsNum(size_t imageID) const
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetNumPoints2D();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
const CKeypoints& CRAMDatabase::GetImageKeypoints(size_t imageID) const
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetKeypoints();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
CKeypoints& CRAMDatabase::GetImageKeypoints(size_t imageID)
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetKeypoints();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
void CRAMDatabase::AddImageKeypoints(size_t imageID, const CKeypoints& keypoints)
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		it->second.SetPoints2D(keypoints);
		return;
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}

size_t CRAMDatabase::GetImageDescriptorsNum(size_t imageID) const
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetNumDescriptors();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
const CSIFTDescriptors& CRAMDatabase::GetImageDescriptors(size_t imageID) const
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetDescriptors();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
CSIFTDescriptors& CRAMDatabase::GetImageDescriptors(size_t imageID)
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		return it->second.GetDescriptors();
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}
void CRAMDatabase::AddImageDescriptors(size_t imageID, const CSIFTDescriptors& descriptors)
{
	auto it = images.find(imageID);
	if (it != images.end())
	{
		it->second.SetDescriptors(descriptors);
		return;
	}
	Check(false, (boost::format("This image does not exist! imageID = %1%") % imageID).str());
}

size_t CRAMDatabase::GetMatchesNum(size_t imageID1, size_t imageID2) const
{
	if (imageID1 > imageID2)
	{
		return GetMatchesNum(imageID2, imageID1);
	}
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 != imageID2);

	auto findFirst = matches.find(imageID1);
	if (findFirst != matches.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second.size();
		}
	}
	return 0;
}
const CSIFTMatches& CRAMDatabase::GetMatches(size_t imageID1, size_t imageID2) const
{
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 < imageID2);

	auto findFirst = matches.find(imageID1);
	if (findFirst != matches.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second;
		}
	}
	Check(false, (boost::format("No matching relationship can be found between Image %1% and Image %1%!") % imageID1 % imageID2).str());
}
CSIFTMatches& CRAMDatabase::GetMatches(size_t imageID1, size_t imageID2)
{
	Check(imageID1 < imageID2);
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());

	auto findFirst = matches.find(imageID1);
	if (findFirst != matches.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second;
		}
	}
	Check(false, (boost::format("No matching relationship can be found between Image %1% and Image %1%!") % imageID1 % imageID2).str());
}
unordered_map<size_t, unordered_map<size_t, CSIFTMatches>> CRAMDatabase::GetAllMatches() const
{
	unordered_map<size_t, unordered_map<size_t, CSIFTMatches>> result;
	result.reserve(matches.size());
	for (const auto& pair : matches)
	{
		unordered_map<size_t, CSIFTMatches> temp;
		temp.reserve(pair.second.size());
		for (const auto& tempPair : pair.second)
		{
			temp[tempPair.first] = tempPair.second;
		}
		result[pair.first] = temp;
	}
	return result;
}
void CRAMDatabase::AddMatches(size_t imageID1, size_t imageID2, const CSIFTMatches& matches)
{
	if (matches.empty())
	{
		return;
	}
	if (imageID1 > imageID2)
	{
		CSIFTMatches matchesReversed = matches;
		for (CSIFTMatch& match : matchesReversed)
		{
			swap(match.point2DIndex1, match.point2DIndex2);
		}
		AddMatches(imageID2, imageID1, matchesReversed);
	}
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 != imageID2);
	if (this->matches[imageID1].find(imageID2) != this->matches[imageID1].end())
	{
		this->matches[imageID1][imageID2].insert(this->matches[imageID1][imageID2].end(), matches.begin(), matches.end());
	}
	else
	{
		this->matches[imageID1][imageID2] = matches;
	}
}
size_t CRAMDatabase::GetMatchedImagesNum(size_t imageID) const
{
	Check(images.find(imageID) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID).str());

	size_t count = 0;
	for (const auto& pair : matches)
	{
		if (pair.first > imageID)
		{
			continue;
		}
		else if (pair.first == imageID)
		{
			count += pair.second.size();
		}
		else if(pair.second.find(imageID) != pair.second.end())
		{
			count++;
		}
	}
	return count;
}

size_t CRAMDatabase::GetInlierMatchesNum(size_t imageID1, size_t imageID2) const
{
	if (imageID1 > imageID2)
	{
		return GetInlierMatchesNum(imageID2, imageID1);
	}
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 != imageID2);

	auto findFirst = twoViewGeometries.find(imageID1);
	if (findFirst != twoViewGeometries.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second.inlierMatches.size();
		}
	}
	return 0;
}
const CTwoViewGeometry& CRAMDatabase::GetTwoViewGeometry(size_t imageID1, size_t imageID2) const
{
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 < imageID2);

	auto findFirst = twoViewGeometries.find(imageID1);
	if (findFirst != twoViewGeometries.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second;
		}
	}
	Check(false, (boost::format("No two-view geometry can be found between Image %1% and Image %1%!") % imageID1 % imageID2).str());
}
CTwoViewGeometry& CRAMDatabase::GetTwoViewGeometry(size_t imageID1, size_t imageID2)
{
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 < imageID2);

	auto findFirst = twoViewGeometries.find(imageID1);
	if (findFirst != twoViewGeometries.end())
	{
		auto findSecond = findFirst->second.find(imageID2);
		if (findSecond != findFirst->second.end())
		{
			return findSecond->second;
		}
	}
	Check(false, (boost::format("No two-view geometry can be found between Image %1% and Image %2%!") % imageID1 % imageID2).str());
}
unordered_map<size_t, unordered_map<size_t, CTwoViewGeometry>> CRAMDatabase::GetAllTwoViewGeometries() const
{
	unordered_map<size_t, unordered_map<size_t, CTwoViewGeometry>> result;
	result.reserve(twoViewGeometries.size());
	for (const auto& pair : twoViewGeometries)
	{
		unordered_map<size_t, CTwoViewGeometry> temp;
		temp.reserve(pair.second.size());
		for (const auto& tempPair : pair.second)
		{
			temp[tempPair.first] = tempPair.second;
		}
		result[pair.first] = temp;
	}
	return result;
}
void CRAMDatabase::AddTwoViewGeometry(size_t imageID1, size_t imageID2, const CTwoViewGeometry& twoViewGeometry)
{
	Check(images.find(imageID1) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID1).str());
	Check(images.find(imageID2) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID2).str());
	Check(imageID1 < imageID2);
	this->twoViewGeometries[imageID1][imageID2] = twoViewGeometry;
	for (const CSIFTMatch& match : twoViewGeometry.inlierMatches)
	{
		images[imageID1].AddCorrespondence(match.point2DIndex1, imageID2, match.point2DIndex2);
		images[imageID2].AddCorrespondence(match.point2DIndex2, imageID1, match.point2DIndex1);
	}
}
size_t CRAMDatabase::GetTwoViewGeometryImagesNum(size_t imageID) const
{
	Check(images.find(imageID) != images.end(), (boost::format("This image does not exist! imageID = %1%") % imageID).str());
	size_t count = 0;
	for (const auto& pair : twoViewGeometries)
	{
		if (pair.first > imageID)
		{
			continue;
		}
		else if (pair.first == imageID)
		{
			count += pair.second.size();
		}
		else if(pair.second.find(imageID) != pair.second.end())
		{
			count++;
		}
	}
	return count;
}
unordered_map<pair<size_t, size_t>, size_t, MatchPairHash, MatchPairEqual> CRAMDatabase::GetAllCorrespondences() const
{
	unordered_map<pair<size_t, size_t>, size_t, MatchPairHash, MatchPairEqual> correspondences;
	correspondences.reserve(twoViewGeometries.size());
	for (const auto& pair : twoViewGeometries)
	{
		const size_t imageID1 = pair.first;
		for (const auto& pair2 : pair.second)
		{
			const size_t imageID2 = pair2.first;
			Check(imageID1 < imageID2);
			correspondences[make_pair(imageID1, imageID2)] = pair2.second.inlierMatches.size();
		}
	}
	return correspondences;
}
void CRAMDatabase::SaveAsDir(const string& dirPath) const
{
	Check(IsDirExists(dirPath));
	const string outputPath = EnsureTrailingSlash(dirPath);

	ofstream cameraFile(outputPath + "cameras.txt");
	Check(cameraFile.is_open());
	cameraFile << "cameras num:" << cameras.size() << endl;
	cameraFile << "cameraID,width,height,cameraModelName,params,isFocalLengthPrior" << endl;
	for (const auto& pair : cameras)
	{
		const CCamera& camera = pair.second;
		const size_t width = camera.GetWidth();
		const size_t height = camera.GetHeight();
		const string& cameraModel = camera.GetCameraModelName();
		const vector<double>& params = camera.GetParams();
		const bool isFocalLengthPrior = camera.IsFocalLengthPrior();

		cameraFile << pair.first << "," << width << "," << height << "," << cameraModel << "," << VectorToString(params, " ") << isFocalLengthPrior << endl;
	}
	cameraFile.close();

	ofstream imageFile(outputPath + "images.txt");
	Check(imageFile.is_open());
	imageFile << "images num:" << images.size() << endl;
	for (const auto& pair : images)
	{


	}






}
