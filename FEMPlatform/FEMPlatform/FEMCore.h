#pragma once
#include "MeshManager/FEMMeshManager.h"
#include "MeshManager/FEM2DMeshManager.h"
#include "MeshManager/FEM3DMeshManager.h"

#include "FEMSolver/FEMSolver.h"
#include "FEMSolver/FEMNRSolveStrategy.h"
#include "FEMSolver/FEMNDDRSolveStrategy.h"
#include "FEMSolver/FEM2DSolver.h"
#include "FEMSolver/FEM2DStaticSolver.h"

#include "MatrixSolver/SluMTMatrixSolver.h"

#include "FEMModel/FEMModel.h"

#include <iostream>
#include <string>
using namespace std;

/*提供用户操作的接口，给出用户在图形界面可以进行的设置*/
class FEMCore
{
public:
	FEMCore();
	~FEMCore();
	void setModel(FEMModel* _model);
	void setAnalysisType(string analysistype);
	void setFEMSolveStrategy(string strategy);
	void setMatrixSolver(string matrixsolver);
	void solve();
	void postprocess();
	void setMaxIterSteps(const int _maxitersteps);
	void setMaxError(const double _error);

private:
	FEMModel* model;
	FEMMeshManager* meshmanager;
	FEMSolver* solver;
};

