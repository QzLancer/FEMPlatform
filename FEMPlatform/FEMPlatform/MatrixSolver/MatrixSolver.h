#pragma once

class MatrixSolver 
{
public:
	virtual ~MatrixSolver() = default;
	virtual double* solveMatrix(int** locs, double* vals, double* F, int valsize, int vecsize) = 0;
};