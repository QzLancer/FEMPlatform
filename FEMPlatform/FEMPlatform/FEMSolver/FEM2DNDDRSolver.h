#pragma once

#include "FEM2DSolver.h"

class FEM2DNDDRSolver :
    public FEM2DSolver
{
public:
    virtual void solveStatic() override;
	void solve2DAxim();
	void solve2DPlane();

private:
	struct NDDR_FEMNode
	{
		int NumberofNeighbourElement{ 0 };
		vector<int> NeighbourElementId;	//和节点相关的单元编号
		vector<int> NeighbourElementNumber;	//节点在对应单元中的编号
	};
	vector<NDDR_FEMNode> nddrnode;
};

