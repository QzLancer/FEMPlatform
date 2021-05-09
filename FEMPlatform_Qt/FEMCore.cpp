#include "FEMCore.h"

FEMCore::FEMCore() :
	model(nullptr),
	meshmanager(nullptr),
	solver(nullptr)
{
	if (solver != nullptr) {
		delete solver;
		solver = nullptr;
	}
	if (meshmanager != nullptr) {
		delete meshmanager;
		meshmanager = nullptr;
	}
}

FEMCore::~FEMCore()
{
}

void FEMCore::setModel(FEMModel* _model)
{
	model = _model;
	//根据维度选择分网管理器类型
	if (model->getDimension() == FEMModel::DIMENSION::D2AXISM || model->getDimension() == FEMModel::DIMENSION::D2PLANE) {
		meshmanager = new FEM2DMeshManager;

	}
	else if (model->getDimension() == FEMModel::DIMENSION::THREE) {
		meshmanager = new FEM3DMeshManager;
	}
	else {
		std::cout << "Error: invalid dimension!\n";
		exit(0);
	}
	if (!model->getGeoFile().empty()) {
		meshmanager->readGeoFile(model->getGeoFile());
		meshmanager->readMeshFile();
	}
	else {
		meshmanager->readMeshFile(model->getMeshFile());
	}
	meshmanager->meshUnitConvert(model->getUnitRatio());
}

void FEMCore::setAnalysisType(string analysistype)
{
	if (analysistype == "static") {
		
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
		solver = new FEM2DNRSolver;
		solver->dimension = model->getDimension();
	}
	else if (solvestrategy == "NDDR") {
		solver = new FEM2DNDDRSolver;
		solver->dimension = model->getDimension();
	}
	else if (solvestrategy == "NDDRGPU") {
		solver = new FEM2DNDDRCUDASolver;
		solver->dimension = model->getDimension();
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

void FEMCore::solveStatic()
{
	solver->setNodes(meshmanager->getNumofNodes(), meshmanager->getNodes());
	solver->setVtxElements(meshmanager->getNumofVtxEle(), meshmanager->getVtxElements());
	solver->setEdgElements(meshmanager->getNumofEdgEle(), meshmanager->getEdgElements());
	solver->setTriElements(meshmanager->getNumofTriEle(), meshmanager->getTriElements());
	solver->setMaterial(model->getMaterialMap());
	solver->setLoad(model->getLoadMap());
	solver->setBoundary(model->getBoundaryMap());
	solver->setDeformedDomain(model->getDeformedList());
	solver->setMovingPart(model->getMovingMap());

	solver->solveStatic();
	//solver->solveMagneticForce();	//电磁力计算，目前是静态特性的思路放在core中调用
}

void FEMCore::postprocess()
{
	solver->writeVtkFile(model->getModelName());
	solver->writeTxtFile(model->getModelName());
}

void FEMCore::setMaxIterSteps(const int _maxitersteps)
{
	solver->setMaxIterSteps(_maxitersteps);
}

void FEMCore::setMaxError(const double _error)
{
	solver->setMaxError(_error);
}
