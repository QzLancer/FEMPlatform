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
		vector<int> NeighbourElementId;	//�ͽڵ���صĵ�Ԫ���
		vector<int> NeighbourElementNumber;	//�ڵ��ڶ�Ӧ��Ԫ�еı��
	};
	vector<NDDR_FEMNode> nddrnode;
};

