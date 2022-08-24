#pragma once

#include "FEMMaterial.h"
#include <vector>
/// @brief 用于存储节点信息，目前整合了NDDR所需的相关信息
struct CNode
{
    /// @brief 节点坐标
    double x{ 0 }, y{ 0 }, z{ 0 };
    /// @brief 边界条件
    int bdr{ 0 };
    /// @brief 磁矢位
    double A{ 0 };
    /// @brief 迭代过程中旧的磁矢位
    double A_old{ 0 };
    /// @brief 二维轴对称等效磁矢位（At要换算到A）
    double At{ 0 };
    /// @brief 二位轴对称迭代过程中旧的等效磁矢位
    double At_old{ 0 };
    /// @brief 二位轴对称非线性迭代的磁矢位变化
    double delta_At{ 0 };
    /// @brief 同上
    double delta_At_old{ 0 };
    /// @brief 节点周围的单元数目
    int NumberofNeighbourElement{ 0 };
    /// @brief 节点周边的单元编号
    int NeighbourElementId[15];
    /// @brief 节点在对应单元中的编号
    int NeighbourElementNumber[15];
    /// @brief 节点电磁力
    double NodeForcex{ 0 }, NodeForcey{ 0 }, NodeForcez{ 0 }, NodeForce{ 0 };
    /// @brief NR迭代过程中计算得到的右侧项
    double SumRHSContri{ 0 };
    /// @brief 周边节点编号
    int NeiborNode[20];
    /// @brief 周边节点数目
    int NumNeiborNodes{ 0 };
    /// @brief 节点周边全部单元对右侧列向量的贡献
    double JsSum{ 0 };
    /// @brief 节点周边全部节点的JsSum之和
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

/// @brief 节点单元
struct CVtxElement {
    int n{0};
    int domain{0};
};

/// @brief 线单元，domain用于边界条件的检索
struct CEdgElement {
    int n[2]{0};
    double x[2]{0};
    double y[2]{0};
    double z[2]{0};
    int domain{0};
};

/// @brief 三角形单元
struct CTriElement {
    int n[3]{0};// ni, nj, nk;//
    double Q[3]{0};// Qi, Qj, Qk;
    double R[3]{0};// Ri, Rj, Rk;
    double C[3][3];// 单元系数矩阵
    double area{0};
    double rc, zc;
    double xdot;
    int domain{0};
    FEMMaterial* material{nullptr};
    double J{0};  //负载，暂时只考虑电流，直流，后续需要单独建一个类
    double Bx{ 0 }, By{ 0 }, Bz{ 0 }, B{ 0 };
    double ElementForcex{ 0 }, ElementForcey{ 0 }, ElementForcez{ 0 }, ElementForce{ 0 };
    double RHSContri[3]{ 0 };
    double mut, mu;
    double ElmRowSum[3][3];
};