#pragma once
#include "FEM2DSolver.h"

//只考虑第一类边界条件计算
class FEM2DStaticSolver :
    public FEM2DSolver
{
    virtual void solve() override;
};

