#pragma once
#include <vector>
#include <tbb/concurrent_vector.h>

#include "../Base/Base.h"


// 一个CTrackElement表示一张影像上的某个2D点对某个3D点的观测
struct CTrackElement
{
	size_t imageID;
	size_t point2DIndex;

	CTrackElement();
	CTrackElement(size_t imageID, size_t point2DIndex);
};


// CTrack表示一个3D点的所有观测. 每个3D点都对应一个CTrack
class CTrack
{
public:
	CTrack() noexcept{};

	inline size_t GetTrackLength() const noexcept
	{
		return elements.size();
	}
	
	inline const std::vector<CTrackElement>& GetAllElements() const
	{
		return elements;
	}
	inline std::vector<CTrackElement>& GetAllElements()
	{
		return elements;
	}
	inline void SetAllElements(const std::vector<CTrackElement>& elements)
	{
		this->elements = elements;
	}
	
	inline const CTrackElement& GetElement(size_t index) const
	{
		Check(index < elements.size());
		return elements[index];
	}
	inline CTrackElement& GetElement(size_t index)
	{
		Check(index < elements.size());
		return elements[index];
	}
	inline void SetElement(size_t index, const CTrackElement& element)
	{
		Check(index < elements.size());
		elements[index] = element;
	}

	inline void AddElement(const CTrackElement& element)
	{
		elements.push_back(element);
	}
	inline void AddElement(size_t imageID, size_t point2DIndex)
	{
		elements.emplace_back(imageID, point2DIndex);
	}
	inline void AddElements(const std::vector<CTrackElement>& elements)
	{
		this->elements.insert(this->elements.end(), elements.begin(), elements.end());
	}

	inline void DeleteElement(size_t index)
	{
		Check(index < elements.size());

		if (index < elements.size() - 1)
		{
			std::move(elements.begin() + index + 1, elements.end(), elements.begin() + index);
		}
		elements.resize(elements.size() - 1);
	}
	void DeleteElement(size_t imageID, size_t point2DIndex);

	inline void Reserve(size_t numElements)
	{
		elements.reserve(numElements);
	}
	inline void Compress()
	{
		elements.shrink_to_fit();
	}

private:
	std::vector<CTrackElement> elements;
};

