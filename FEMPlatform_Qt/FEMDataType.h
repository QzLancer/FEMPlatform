#pragma once

#include "FEMMaterial.h"
#include <vector>
//Ŀǰ������NDDRNODE��ص����ݣ�������Ҫ�޸�
struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
    int bdr{ 0 }; //�߽�����
    double A{ 0 };
    double A_old{ 0 };
    double At{ 0 };
    double At_old{ 0 };
    int NumberofNeighbourElement{ 0 };
    int NeighbourElementId[15];	//�ͽڵ���صĵ�Ԫ���
    int NeighbourElementNumber[15];	//�ڵ��ڶ�Ӧ��Ԫ�еı��
    double NodeForcex{ 0 }, NodeForcey{ 0 }, NodeForcez{ 0 }, NodeForce{ 0 };
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
    double C[3][3];// ��Ԫϵ������
    double area{0};
    double rc, zc;
    double xdot;
    int domain{0};
    FEMMaterial* material{nullptr};
    double J{0};  //���أ���ʱֻ���ǵ�����ֱ����������Ҫ������һ����
    double Bx{ 0 }, By{ 0 }, Bz{ 0 }, B{ 0 };
    double ElementForcex{ 0 }, ElementForcey{ 0 }, ElementForcez{ 0 }, ElementForce{ 0 };
    double RHSContri{ 0 };
};