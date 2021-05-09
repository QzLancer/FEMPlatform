#pragma once

#include <vector>

using namespace std;

class MatrixSolver 
{
public:
	virtual ~MatrixSolver() = default;
	virtual vector<double> solveMatrix(vector<vector<int>> locs, vector<double> vals, vector<double> F, int valsize, int vecsize) = 0;
};