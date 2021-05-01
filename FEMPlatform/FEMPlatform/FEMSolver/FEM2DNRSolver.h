#pragma once

#include "FEM2DSolver.h"

#include <vector>
#include <math.h>

//只考虑第一类边界条件计算
class FEM2DNRSolver :
    public FEM2DSolver
{
public:
    virtual void solveStatic() override;
    void solve2DAxim();
    void solve2DPlane();
protected:


private:
    
};

