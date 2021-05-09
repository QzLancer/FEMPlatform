#pragma once
#include "MeshManager/FEMMeshManager.h"
#include "MeshManager/FEM2DMeshManager.h"
#include "MeshManager/FEM3DMeshManager.h"

#include "FEMSolver/FEMSolver.h"
#include "FEMSolver/FEM2DSolver.h"
#include "FEMSolver/FEM2DNRSolver.h"
#include "FEMSolver/FEM2DNDDRSolver.h"
#include "FEMSolver/FEM2DNDDRCUDASolver.cuh"

#include "MatrixSolver/SluMTMatrixSolver.h"

#include "FEMModel/FEMModel.h"

#include <iostream>
#include <string>
using namespace std;

/*�ṩ�û������Ľӿڣ������û���ͼ�ν�����Խ��е�����*/
class FEMCore
{
public:
	FEMCore();
	~FEMCore();
	void setModel(FEMModel* _model);
	void setAnalysisType(string analysistype);
	void setFEMSolveStrategy(string strategy);
	void setMatrixSolver(string matrixsolver);
	void solveStatic();
	void postprocess();
	void setMaxIterSteps(const int _maxitersteps);
	void setMaxError(const double _error);

private:
	FEMModel* model;
	FEMMeshManager* meshmanager;
	FEMSolver* solver;
};

