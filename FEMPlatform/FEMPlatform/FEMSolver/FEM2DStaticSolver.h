#pragma once

#include "FEM2DSolver.h"

#include <vector>
#include <algorithm>

//只考虑第一类边界条件计算
class FEM2DStaticSolver :
    public FEM2DSolver
{
public:
    virtual void solve() override;

private:
    
};

