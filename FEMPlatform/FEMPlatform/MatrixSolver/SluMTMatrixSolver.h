#pragma once
#include "MatrixSolver.h"

#include <slu_mt_ddefs.h>
#include <map>
#include <thread>

class SluMTMatrixSolver :
    public MatrixSolver
{
public:
    SluMTMatrixSolver();
    virtual ~SluMTMatrixSolver() override;
    virtual double* solveMatrix(vector<vector<int>> locs, vector<double> vals, vector<double> F, int valsize, int vecsize) override;
    void setNumberofProcs(const int _nprocs);

private:
    SuperMatrix A, B, L, U;
    int nprocs;
};

