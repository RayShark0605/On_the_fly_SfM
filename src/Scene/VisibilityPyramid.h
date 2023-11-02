#pragma once
#include <vector>
#include <Eigen/Core>

#include "../Base/Base.h"
#include "../Base/Math.h"

// 捕捉2D网格中可见的3D点分布的类. 通过一个得分来捕捉点的分布情况, 更高的得分表示点在网格中的分布更均匀.
// 这个得分是通过多分辨率金字塔中填充了元素的单元格数量来计算的. 如果一个单元格包含至少一个点, 它就会对总得分做出贡献, 而贡献的分数是根据它在金字塔中的分辨率来确定的, 在更高分辨率层级的单元格将对总得分作出更高的贡献。
class CVisibilityPyramid final
{
public:
	CVisibilityPyramid() noexcept
	{
		width = 0;
		height = 0;
		score = 0;
		maxScore = 0;
		pyramid.clear();
	}
	CVisibilityPyramid(size_t numLevels, size_t width, size_t height) noexcept
	{
		this->width = width;
		this->height = height;
		score = 0;
		maxScore = 0;
		pyramid.resize(numLevels);
		for (size_t level = 0; level < numLevels; level++)
		{
			const size_t dim = 1 << (level + 1);
			pyramid[level].setZero(dim, dim);
			maxScore += dim * dim * dim * dim;
		}
	}
	void SetPoint(double x, double y)
	{
		Check(pyramid.size() > 0);
		size_t cx = 0;
		size_t cy = 0;
		CellForPoint(x, y, &cx, &cy);
		for (int i = pyramid.size() - 1; i >= 0; i--)
		{
			Eigen::MatrixXi& level = pyramid[i];
			level(cy, cx) += 1;
			if (level(cy, cx) == 1) 
			{
				score += level.size();
			}
			cx = cx >> 1;
			cy = cy >> 1;
		}
		Check(score <= maxScore);
	}
	void ResetPoint(double x, double y)
	{
		Check(pyramid.size() > 0);
		size_t cx = 0;
		size_t cy = 0;
		CellForPoint(x, y, &cx, &cy);
		for (int i = pyramid.size() - 1; i >= 0; i--)
		{
			Eigen::MatrixXi& level = pyramid[i];
			level(cy, cx) -= 1;
			if (level(cy, cx) == 0) 
			{
				score -= level.size();
			}
			cx = cx >> 1;
			cy = cy >> 1;
		}
		Check(score <= maxScore);
	}

	inline size_t NumLevels() const noexcept
	{
		return pyramid.size();
	}
	inline size_t Width() const noexcept
	{
		return width;
	}
	inline size_t Height() const noexcept
	{
		return height;
	}
	inline size_t Score() const noexcept
	{
		return score;
	}
	inline size_t MaxScore() const noexcept
	{
		return maxScore;
	}

	inline std::string WriteToString() const
	{
		std::ostringstream oss;
		oss << width << "," << height << "," << score << "," << maxScore << "," << pyramid.size() << ",";
		for (const auto& matrix : pyramid) 
		{
			oss << "{";
			for (int i = 0; i < matrix.rows(); i++) 
			{
				for (int j = 0; j < matrix.cols(); j++) 
				{
					oss << matrix(i, j);
					if (j < matrix.cols() - 1) oss << ",";
				}
				if (i < matrix.rows() - 1) oss << ";";
			}
			oss << "}";
		}
		return oss.str();
	}
	inline void ReadFromString(const std::string& str)
	{
		std::istringstream iss(str);
		char ch;
		size_t numLevels;
		iss >> width >> ch >> height >> ch >> score >> ch >> maxScore >> ch >> numLevels >> ch;
		pyramid.resize(numLevels);
		for (size_t level = 0; level < numLevels; level++)
		{
			iss >> ch;  // 读取 '{'
			std::vector<std::vector<int>> values;
			std::vector<int> row;
			int num;
			char sep;

			while (iss >> num)
			{
				row.push_back(num);
				iss >> sep;
				if (sep == ';')
				{
					values.push_back(row);
					row.clear();
				}
				else if (sep == '}')
				{
					values.push_back(row);
					break;
				}
				else if (sep != ',')
				{
					// 错误处理
					return;
				}
			}

			Eigen::MatrixXi matrix(values.size(), values[0].size());
			for (int i = 0; i < values.size(); ++i)
			{
				for (int j = 0; j < values[0].size(); ++j)
				{
					matrix(i, j) = values[i][j];
				}
			}
			pyramid[level] = matrix;
		}
	}

private:
	size_t width, height;
	size_t score;    // 当前的可见性得分
	size_t maxScore; // 当所有单元格都填充时的最大得分
	std::vector<Eigen::MatrixXi> pyramid; // 存储可见性金字塔. 每一层是一个整型矩阵

	void CellForPoint(double x, double y, size_t* cx, size_t* cy) const noexcept
	{
		Check(width > 0 && height > 0);
		const size_t maxDim = 1 << pyramid.size();
		*cx = Clamp<size_t>(maxDim * x / width, 0, maxDim - 1);
		*cy = Clamp<size_t>(maxDim * y / height, 0, maxDim - 1);
	}

};








