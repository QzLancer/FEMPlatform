#pragma once

#include "FEM2DSolver.h"

class FEM2DNDDRSolver :
    public FEM2DSolver
{
public:
    virtual void solveStatic() override;
	void solve2DAxim();
	void solve2DAxim1();
	void solve2DAxim2();
	void solve2DPlane();
	void solve2DPlane1();
	void solve2DPlane2();
	void Update_Magnetic_Node_A();
	void Update_Magnetic_Node_A_old();
private:
	struct NDDR_FEMNode
	{
		int NumberofNeighbourElement{ 0 };
		vector<int> NeighbourElementId;	//�ͽڵ���صĵ�Ԫ���
		vector<int> NeighbourElementNumber;	//�ڵ��ڶ�Ӧ��Ԫ�еı��
	};
};

