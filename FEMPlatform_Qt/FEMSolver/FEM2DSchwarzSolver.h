#pragma once
#include "FEM2DSolver.h"
#include "mpmetis.h"

#include <string>

class FEM2DSchwarzSolver :
    public FEM2DSolver
{
public:
    FEM2DSchwarzSolver(int _numofdomain = 0);
    ~FEM2DSchwarzSolver();

    virtual void solveStatic() override;
    virtual void solveDynamic() override;

protected:
    virtual void processBoundaryCondition() override;

private:
    void generateMetisMesh();
    void solve2DAximLinear();
    void solve2DAximNonlinear();

    char metisfile[256];
    int numofdomain;    //��������
    int* d_nodesize;  //�����ڵĽڵ���Ŀ
    int* d_dof;   //�����ڵ����ɽڵ���Ŀ
    int* d_elesize;  //�����ڵĵ�Ԫ����
    int** d_eleid;  //�ֲ���Ԫ���ӳ�䵽ȫ�ֵ�Ԫ���
    int** eparttable;   //ȫ�ֵ�Ԫ��ŵ��ֲ���Ԫ���
    int** crossbdrtable; //���紦�߽�Ľڵ���
    int** nparttable;   //ȫ�ֽڵ���ӳ�䵽�ֲ��ڵ���
    int** d_nodeid; //�ֲ��ڵ���ӳ�䵽ȫ�ֽڵ���
    int** d_node_reorder; //�Ծֲ��ڵ㰴�����ɽڵ�-�߽�ڵ��˳���������򣬽��������ڵ��ŵ�δ����ڵ��ŵ�ӳ��
    int** d_node_pos;   //����δ����ڵ��ŵ������ڵ��ŵ�ӳ��
};

