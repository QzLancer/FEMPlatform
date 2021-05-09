#pragma once
#include "FEMSolver.h"
class FEM2DSolver :
    public FEMSolver
{
public:
    virtual void solveStatic() = 0;
    virtual void solveDynamic() = 0;
    virtual void solveMagneticForce() override;

protected:
    void makeTrangle();
    virtual void processBoundaryCondition() override;    //处理自由节点
    virtual void processMaterial() override;    //将材料添加到单元上
    virtual void processLoad() override;    //将电流负载添加到节点上
    void updateB();
    void updateB(int i_tri);
};

