#pragma once
struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
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
    double P[3]{0};// Pi, Pj, Pk;
    double Q[3]{0};// Qi, Qj, Qk;
    double area;
    double rc, zc;
    double ydot;
    int domain{0};
};