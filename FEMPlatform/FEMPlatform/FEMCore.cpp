#include "FEMCore.h"

void FEMCore::setDimension(string dimension)
{
	if (dimension == "2D") {
		generator = new FEM2DGenerator;
		
	}
	else if (dimension == "3D") {
		generator = new FEM3DGenerator;
	}
	else {
		std::cout << "Error: invalid dimension!";
		exit(0);
	}
}

void FEMCore::readMeshData(string meshfile)
{
	generator->readMeshFile(meshfile);
}

void FEMCore::createElement2Material()
{
	generator->createElement2Material();
}

void FEMCore::bulidGeometry2Load()
{
	generator->bulidGeometry2Load();
}

void FEMCore::buildGeometry2Constrain()
{
	generator->buildGeometry2Constrain();
}

void FEMCore::setAnalysisType(string analysistype)
{
	if (analysistype == "static") {
		solver = new FEM2DStaticSolver;
	}
	else if (analysistype == "dynamic") {
		
	}
}

void FEMCore::setFEMSolver(string solvertype)
{
	if (solvertype == "NR") {
		solver->setSolveStrategy(new FEMNRSolveStrategy);
	}
	else if (solvertype == "NDDR") {
		solver->setSolveStrategy(new FEMNDDRSolveStrategy);
	}
}

void FEMCore::solve()
{
	//solver.setdata(generator.getdata());
}

void FEMCore::postoperation()
{
}
