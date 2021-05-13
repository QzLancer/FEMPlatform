#pragma once
#include "FEMSolver.h"
class FEM2DSolver :
    public FEMSolver
{
public:
    virtual void solveStatic() = 0;
    virtual void solveDynamic() = 0;
    virtual void solveMagneticForce() override;
    virtual void solveMagneticForce1() override;
    virtual double solveSpringForce(int domain, double pos);
    virtual double solveEnergy();
    virtual double solveInductance(double I);
    virtual double solveFlux(int domain);

protected:
    void makeTrangle();
    virtual void processBoundaryCondition() override;    //处理自由节点
    virtual void processMaterial() override;    //将材料添加到单元上
    virtual void processLoad() override;    //将电流负载添加到节点上
    void setCoilCurrent(double I);
    void updateB();
    void updateB(int i_tri);
    double Fx, Fy;
};

