#pragma once
#include "FEMSolver.h"
class FEM2DSolver :
    public FEMSolver
{
public:
    virtual void solve() = 0;

protected:
    void makeTrangle();
    virtual void processBoundaryCondition() override;    //�������ɽڵ�
    virtual void processMaterial() override;    //��������ӵ���Ԫ��
    virtual void processLoad() override;    //������������ӵ��ڵ���
    void updateB();
};

