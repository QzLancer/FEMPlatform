#pragma once
#include "FEMGenerator.h"
#include "FEM2DGenerator.h"
#include "FEM3DGenerator.h"
#include "FEMSolver/FEMSolver.h"
#include "FEMSolver/FEMNRSolveStrategy.h"
#include "FEMSolver/FEMNDDRSolveStrategy.h"
#include "FEMSolver/FEM2DSolver.h"
#include "FEMSolver/FEM2DStaticSolver.h"
#include "MatrixSolver/SluMTMatrixSolver.h"

#include <iostream>
#include <string>
using namespace std;

/*提供用户操作的接口，给出用户在图形界面可以进行的设置*/
class FEMCore
{
public:
	FEMCore();
	~FEMCore();
	void setDimension(string dimension);
	void readMeshData(string meshfile);
	void createElement2Material();
	void bulidGeometry2Load();
	void buildGeometry2Constrain();
	void setAnalysisType(string analysistype);
	void setFEMSolveStrategy(string strategy);
	void setMatrixSolver(string matrixsolver);
	void solve();
	void postoperation();
	void setMaxIterSteps(const int _maxitersteps);

private:
	FEMGenerator* generator;
	FEMSolver* solver;

};

