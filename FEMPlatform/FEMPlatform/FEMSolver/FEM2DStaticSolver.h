#pragma once
#include "FEM2DSolver.h"

//ֻ���ǵ�һ��߽���������
class FEM2DStaticSolver :
    public FEM2DSolver
{
    virtual void solve() override;
};

