#pragma once
#include "FEMSolver.h"
class FEM2DSolver :
    public FEMSolver
{
public:
    virtual void solve() = 0;

protected:
    void makeTrangle();
    virtual void processBoundaryCondition() override;    //处理自由节点
    virtual void processMaterial() override;    //将材料添加到单元上
    virtual void processLoad() override;    //将电流负载添加到节点上
    void updateB();
};

