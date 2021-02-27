#include "FEMCore.h"

FEMCore::FEMCore() :
	generator(nullptr),
	solver(nullptr)
{
	if (solver != nullptr)
		delete solver;
	if (generator != nullptr)
		delete generator;
}

FEMCore::~FEMCore()
{
}

void FEMCore::setDimension(string dimension)
{
	if (dimension == "2D") {
		generator = new FEM2DGenerator;
		
	}
	else if (dimension == "3D") {
		generator = new FEM3DGenerator;
	}
	else {
		std::cout << "Error: invalid dimension!\n";
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
	else {
		cout << "Error: invalid analysis type!\n";
		exit(0);
	}
}

void FEMCore::setFEMSolveStrategy(string solvestrategy)
{
	if (solvestrategy == "NR") {
		solver->setSolveStrategy(new FEMNRSolveStrategy);
	}
	else if (solvestrategy == "NDDR") {
		solver->setSolveStrategy(new FEMNDDRSolveStrategy);
	}
	else {
		cout << "Error: invalid solve strategy!\n";
		exit(0);
	}
}

void FEMCore::setMatrixSolver(string matrixsolver)
{
	if (matrixsolver == "SuperLU_MT") {
		solver->setMatrixSolver(new SluMTMatrixSolver);
	}
	else {
		cout << "Error: invalid matrix solver!\n";
		exit(0);
	}
}

void FEMCore::solve()
{
	solver->setNodes(generator->getNumofNodes(), generator->getNodes());
	solver->setVtxElements(generator->getNumofVtxEle(), generator->getVtxElements());
	solver->setEdgElements(generator->getNumofEdgEle(), generator->getEdgElements());
	solver->setTriElements(generator->getNumofTriEle(), generator->getTriElements());
	solver->setMaterial(generator->getMaterial());
	solver->setLoad(generator->getLoad());
	solver->setBoundary(generator->getBoundary());
	solver->solve();
}

void FEMCore::postoperation()
{
}

void FEMCore::setMaxIterSteps(const int _maxitersteps)
{
	solver->setMaxIterSteps(_maxitersteps);
}
