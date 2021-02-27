#pragma once
#include "FEMSolver.h"
class FEM2DSolver :
    public FEMSolver
{
public:
    virtual void solve() = 0;

protected:
    void makeTrangle(int index);
};

