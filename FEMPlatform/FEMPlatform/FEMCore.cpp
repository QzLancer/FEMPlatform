#include "FEMCore.h"

FEMCore::FEMCore() :
	model(nullptr),
	meshmanager(nullptr),
	solver(nullptr)
{
	if (solver != nullptr)
		delete solver;
	if (meshmanager != nullptr)
		delete meshmanager;
}

FEMCore::~FEMCore()
{
}

void FEMCore::setModel(FEMModel* _model)
{
	model = _model;
	//根据维度选择分网管理器类型
	if (model->getDimension() == FEMModel::DIMENSION::TWO) {
		meshmanager = new FEM2DMeshManager;

	}
	else if (model->getDimension() == FEMModel::DIMENSION::THREE) {
		meshmanager = new FEM3DMeshManager;
	}
	else {
		std::cout << "Error: invalid dimension!\n";
		exit(0);
	}
	meshmanager->readMeshFile(model->getMeshFile());
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
	solver->setNodes(meshmanager->getNumofNodes(), meshmanager->getNodes());
	solver->setVtxElements(meshmanager->getNumofVtxEle(), meshmanager->getVtxElements());
	solver->setEdgElements(meshmanager->getNumofEdgEle(), meshmanager->getEdgElements());
	solver->setTriElements(meshmanager->getNumofTriEle(), meshmanager->getTriElements());
	solver->setMaterial(model->getMaterialMap());
	solver->setLoad(model->getLoadMap());
	solver->setBoundary(model->getBoundaryMap());
	solver->solve();
}

void FEMCore::postprocess()
{
	solver->writeVtkFile(model->getModelName());
}

void FEMCore::setMaxIterSteps(const int _maxitersteps)
{
	solver->setMaxIterSteps(_maxitersteps);
}

void FEMCore::setMaxError(const double _error)
{
	solver->setMaxError(_error);
}
