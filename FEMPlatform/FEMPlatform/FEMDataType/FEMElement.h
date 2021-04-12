#pragma once

#include "../FEMMaterial.cuh"

struct CVtxElement {
    int n{ 0 };
    int domain{ 0 };
};

struct CEdgElement {
    int n[2]{ 0 };
    double x[2]{ 0 };
    double y[2]{ 0 };
    double z[2]{ 0 };
    int domain{ 0 };
};

class CTriElement {
public:
    void makeTrangle();

    int n[3]{ 0 };// ni, nj, nk;//
    double Q[3]{ 0 };// Qi, Qj, Qk;
    double R[3]{ 0 };// Ri, Rj, Rk;
    double C[3][3];// ��Ԫϵ������
    double area{ 0 };
    double rc, zc;
    double xdot;
    int domain{ 0 };
    FEMMaterial* material{ new FEMMaterial };
    double J{ 0 };  //���أ���ʱֻ���ǵ�����ֱ����������Ҫ������һ����
};