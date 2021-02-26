#pragma once
struct C2DNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
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
    double x[3]{0};
    double y[3]{0};
    double z[3]{0};
    int domain{0};
};