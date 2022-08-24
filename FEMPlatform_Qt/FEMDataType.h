#pragma once

#include "FEMMaterial.h"
#include <vector>
/// @brief ���ڴ洢�ڵ���Ϣ��Ŀǰ������NDDR����������Ϣ
struct CNode
{
    /// @brief �ڵ�����
    double x{ 0 }, y{ 0 }, z{ 0 };
    /// @brief �߽�����
    int bdr{ 0 };
    /// @brief ��ʸλ
    double A{ 0 };
    /// @brief ���������оɵĴ�ʸλ
    double A_old{ 0 };
    /// @brief ��ά��ԳƵ�Ч��ʸλ��AtҪ���㵽A��
    double At{ 0 };
    /// @brief ��λ��ԳƵ��������оɵĵ�Ч��ʸλ
    double At_old{ 0 };
    /// @brief ��λ��ԳƷ����Ե����Ĵ�ʸλ�仯
    double delta_At{ 0 };
    /// @brief ͬ��
    double delta_At_old{ 0 };
    /// @brief �ڵ���Χ�ĵ�Ԫ��Ŀ
    int NumberofNeighbourElement{ 0 };
    /// @brief �ڵ��ܱߵĵ�Ԫ���
    int NeighbourElementId[15];
    /// @brief �ڵ��ڶ�Ӧ��Ԫ�еı��
    int NeighbourElementNumber[15];
    /// @brief �ڵ�����
    double NodeForcex{ 0 }, NodeForcey{ 0 }, NodeForcez{ 0 }, NodeForce{ 0 };
    /// @brief NR���������м���õ����Ҳ���
    double SumRHSContri{ 0 };
    /// @brief �ܱ߽ڵ���
    int NeiborNode[20];
    /// @brief �ܱ߽ڵ���Ŀ
    int NumNeiborNodes{ 0 };
    /// @brief �ڵ��ܱ�ȫ����Ԫ���Ҳ��������Ĺ���
    double JsSum{ 0 };
    /// @brief �ڵ��ܱ�ȫ���ڵ��JsSum֮��
    double SumNeiborJsSum{ 0 };
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

/// @brief �ڵ㵥Ԫ
struct CVtxElement {
    int n{0};
    int domain{0};
};

/// @brief �ߵ�Ԫ��domain���ڱ߽������ļ���
struct CEdgElement {
    int n[2]{0};
    double x[2]{0};
    double y[2]{0};
    double z[2]{0};
    int domain{0};
};

/// @brief �����ε�Ԫ
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
    double RHSContri[3]{ 0 };
    double mut, mu;
    double ElmRowSum[3][3];
};