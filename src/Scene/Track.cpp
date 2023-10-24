#include "Track.h"
#include <limits>

using namespace std;

CTrackElement::CTrackElement()
{
	imageID = numeric_limits<size_t>::max();
	point2DIndex = numeric_limits<size_t>::max();
}
CTrackElement::CTrackElement(size_t imageID, size_t point2DIndex)
{
	this->imageID = imageID;
	this->point2DIndex = point2DIndex;
}
void CTrack::DeleteElement(size_t imageID, size_t point2DIndex)
{
	elements.erase(remove_if(elements.begin(), elements.end(), [imageID, point2DIndex](const CTrackElement& element)
		{
			return element.imageID == imageID && element.point2DIndex == point2DIndex;
		}), elements.end());
}
