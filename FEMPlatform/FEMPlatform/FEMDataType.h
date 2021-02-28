#pragma once

#include "FEMMaterial.h"

struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
    int bdr{ 0 }; //边界条件
    double A{ 0 };
};

struct C2DNode:
    public CNode
{
    //double y{ 0 };
};

struct C3DNode:
    public C2DNode
{
    //double z{ 0 };
};

struct CVtxElement {
    int n{0};
    int domain{0};
};

struct CEdgElement {
    int n[2]{0};
    double x[2]{0};
    double y[2]{0};
    double z[2]{0};
    int domain{0};
};

struct CTriElement {
    int n[3]{0};// ni, nj, nk;//
    double Q[3]{0};// Qi, Qj, Qk;
    double R[3]{0};// Ri, Rj, Rk;
    double C[3][3];// 单元系数矩阵
    double area{0};
    double rc, zc;
    double ydot;
    int domain{0};
    FEMMaterial* material{new FEMMaterial};
    double J{0};  //负载，暂时只考虑电流，直流，后续需要单独建一个类
};