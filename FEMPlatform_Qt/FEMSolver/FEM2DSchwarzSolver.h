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
    int numofdomain;    //子域数量
    int* d_nodesize;  //子域内的节点数目
    int* d_dof;   //子域内的自由节点数目
    int* d_elesize;  //子域内的单元数量
    int** d_eleid;  //局部单元编号映射到全局单元编号
    int** eparttable;   //全局单元编号到局部单元编号
    int** crossbdrtable; //交界处边界的节点编号
    int** nparttable;   //全局节点编号映射到局部节点编号
    int** d_nodeid; //局部节点编号映射到全局节点编号
    int** d_node_reorder; //对局部节点按照自由节点-边界节点的顺序重新排序，建立排序后节点编号到未排序节点编号的映射
    int** d_node_pos;   //建立未排序节点编号到排序后节点编号的映射
};

