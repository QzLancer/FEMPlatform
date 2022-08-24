#pragma once

#include <vector>

using namespace std;

class MatrixSolver 
{
public:
	virtual ~MatrixSolver() = default;

	/// @brief [S][A]=[F]�������Ľӿ�
	/// @param locs [S]�з�0Ԫ�ص��������������
	/// @param vals [S]�з�0Ԫ�ص�ֵ
	/// @param F [F]��ֵ
	/// @param valsize [S]�з�0Ԫ�ص���Ŀ��ͬһλ�ÿ��ܳ��ֶ��Ԫ�أ�
	/// @param vecsize [F]�ĳ��� 
	/// @return [A]��ֵ
	virtual vector<double> solveMatrix(vector<vector<int>> locs, vector<double> vals, vector<double> F, int valsize, int vecsize) = 0;
};