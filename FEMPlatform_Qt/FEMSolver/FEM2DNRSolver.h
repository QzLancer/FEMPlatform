#pragma once

#include "FEM2DSolver.h"

#include <vector>
#include <math.h>

//ֻ���ǵ�һ��߽���������
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

