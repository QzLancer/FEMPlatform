#pragma once
#include "FEMGenerator.h"
#include "FEM2DGenerator.h"
#include "FEM3DGenerator.h"
#include "FEMSolver/FEMSolver.h"
#include "FEMSolver/FEMNRSolveStrategy.h"
#include "FEMSolver/FEMNDDRSolveStrategy.h"
#include "FEMSolver/FEM2DSolver.h"
#include "FEMSolver/FEM2DStaticSolver.h"

#include <iostream>
#include <string>
using namespace std;

/*�ṩ�û������Ľӿڣ������û���ͼ�ν�����Խ��е�����*/
class FEMCore
{
public:
	void setDimension(string dimension);
	void readMeshData(string meshfile);
	void createElement2Material();
	void bulidGeometry2Load();
	void buildGeometry2Constrain();
	void setAnalysisType(string analysistype);
	void setFEMSolver(string solvertype);
	void solve();
	void postoperation();

private:
	FEMGenerator* generator;
	FEMSolver* solver;

};

