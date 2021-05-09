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
    virtual void processBoundaryCondition() override;    //�������ɽڵ�
    virtual void processMaterial() override;    //���������ӵ���Ԫ��
    virtual void processLoad() override;    //�������������ӵ��ڵ���
    void updateB();
    void updateB(int i_tri);
};
