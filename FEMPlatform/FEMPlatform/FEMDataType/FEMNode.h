#pragma once
struct CNode
{
    double x{ 0 }, y{ 0 }, z{ 0 };
    int bdr{ 0 }; //±ß½çÌõ¼þ
    double A{ 0 };
};

struct C2DNode :
    public CNode
{
    //double y{ 0 };
};

struct C3DNode :
    public C2DNode
{
    //double z{ 0 };
};