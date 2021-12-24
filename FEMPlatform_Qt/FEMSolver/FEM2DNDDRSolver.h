#pragma once

#include "FEM2DSolver.h"

class FEM2DNDDRSolver :
    public FEM2DSolver
{
public:
    virtual void solveStatic() override;
	virtual void solveDynamic() override;
	void solve2DAxim();
	void solve2DAxim1();
	void solve2DAxim2();
	void solve2DAximOpt();
	void solve2DAximOpt1();
	void solve2DAximPrecondition();
	void solve2DAximRobin();
	void solve2DPlane();
	void solve2DPlane1();
	void solve2DPlane2();
	void Update_Magnetic_Node_A();
	void Update_Magnetic_Node_A_old();
private:
	struct NDDR_FEMNode
	{
		int NumberofNeighbourElement{ 0 };
		vector<int> NeighbourElementId;	//和节点相关的单元编号
		vector<int> NeighbourElementNumber;	//节点在对应单元中的编号
	};

	void DataPrepare();
	void JsSumCalculate();
	void SumNeiborJsSumCalculate();
	void ElmRHSContriCalculate();
	void SumNodeRHSCalculate();
	void UpdateSolutiontoA1();
	void CopyA1toA0();
	void UpdateVe();

	double Gamma = 5;
};

