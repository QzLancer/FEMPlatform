#pragma once

#include <vector>

using namespace std;

class MatrixSolver 
{
public:
	virtual ~MatrixSolver() = default;

	/// @brief [S][A]=[F]矩阵求解的接口
	/// @param locs [S]中非0元素的行坐标和列坐标
	/// @param vals [S]中非0元素的值
	/// @param F [F]的值
	/// @param valsize [S]中非0元素的数目（同一位置可能出现多个元素）
	/// @param vecsize [F]的长度 
	/// @return [A]的值
	virtual vector<double> solveMatrix(vector<vector<int>> locs, vector<double> vals, vector<double> F, int valsize, int vecsize) = 0;
};