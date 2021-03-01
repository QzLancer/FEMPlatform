#include "FEMSolver.h"

FEMSolver::FEMSolver():
    strategy(nullptr),
    matsolver(nullptr),
    m_num_nodes(0),
    m_num_vtxele(0),
    m_num_edgele(0),
    m_num_triele(0),
    mp_node(nullptr),
    mp_vtxele(nullptr),
    mp_edgele(nullptr),
    mp_triele(nullptr),
    maxitersteps(100)
{

}

FEMSolver::~FEMSolver()
{
    if (mp_triele != nullptr)
        delete[] mp_triele;
    if (mp_edgele != nullptr)
        delete[] mp_edgele;
    if (mp_vtxele != nullptr)
        delete[] mp_vtxele;
    if (mp_node != nullptr)
        delete[] mp_node;
    if(strategy != nullptr)
        delete strategy;
}

void FEMSolver::setSolveStrategy(FEMSolveStrategy* _strategy)
{
    strategy = _strategy;
}

void FEMSolver::setMatrixSolver(MatrixSolver* const _matsolver)
{
    matsolver = _matsolver;
}

void FEMSolver::setNodes(const int _numofnodes, CNode* const _nodes)
{
    m_num_nodes = _numofnodes;
    mp_node = _nodes;
}

void FEMSolver::setVtxElements(const int _numofvtx, CVtxElement* const _vtxele)
{
    m_num_vtxele = _numofvtx;
    mp_vtxele = _vtxele;
}

void FEMSolver::setEdgElements(const int _numofedg, CEdgElement* const _edgele)
{
    m_num_edgele = _numofedg;
    mp_edgele = _edgele;
}

void FEMSolver::setTriElements(const int _numoftri, CTriElement* const _triele)
{
    m_num_triele = _numoftri;
    mp_triele = _triele;
}

void FEMSolver::setMaterial(const std::map<int, FEMMaterial*> _materialmap)
{
    materialmap = _materialmap;
}

void FEMSolver::setLoad(const std::map<int, double> _loadmap)
{
    loadmap = _loadmap;
}

void FEMSolver::setBoundary(const std::map<int, FEMBoundary*> _boundarymap)
{
    boundarymap = _boundarymap;
}

void FEMSolver::setMaxIterSteps(const int _maxitersteps)
{
    maxitersteps = _maxitersteps;
}

void FEMSolver::writeVtkFile(std::string _name)
{
	std::string name = std::string("../../result/") + _name + std::string(".vtk");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: openning file!" << endl;
		exit(0);
	}
	/*
		 1: points
		 3: line
		 5: Triangular element
		 9: Quadrilateral element
		10: Tetrahedral element
		12: Hexahedral element
		13: Triangular prism element
		14: Pyramid element
	*/
	/** 数据版本声明 **/
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	/** 标题 **/
	fprintf(fp, "vtk title\n");
	/** 文件格式声明 **/
	fprintf(fp, "ASCII\n");
	/** 几何拓扑结构 **/
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	//节点
	fprintf(fp, "\nPOINTS %d float\n", m_num_nodes);
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%lf %lf %lf\n", mp_node[i].x, mp_node[i].y, mp_node[i].z);
	}
	//一阶三角形单元
	fprintf(fp, "\nCELLS %d %d\n", m_num_triele, 4 * m_num_triele);
	for (int i = 0; i < m_num_triele; ++i) {
		fprintf(fp, "3 %d %d %d\n", mp_triele[i].n[0], mp_triele[i].n[1], mp_triele[i].n[2]);
	}
	fprintf(fp, "\nCELL_TYPES %d\n", m_num_triele);
	int type = 5;
	for (int i = 0; i < m_num_triele; ++i) {
		fprintf(fp, "%d\n", type);
	}
	//节点磁矢位
	fprintf(fp, "\nPOINT_DATA %d\n", m_num_nodes);
	fprintf(fp, "SCALARS A double 1\n");
	fprintf(fp, "LOOKUP_TABLE %s\n", "Atable");
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%f\n", A[i]);
	}

	//单元标量磁感应强度
	fprintf(fp, "\nCELL_DATA %d\n", m_num_triele);
	fprintf(fp, "SCALARS %s double %d\n", "Bnorm", 1);
	fprintf(fp, "LOOKUP_TABLE %s\n", "Btable");
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		fprintf(fp, "%lf\n", B[i_tri]);
	}

	fprintf(fp, "\nVECTORS %s double\n", "Bvector");
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		fprintf(fp, "%lf %lf %lf\n", Bx[i_tri], By[i_tri], 0.0);
	}
	fclose(fp);
}

std::vector<double> FEMSolver::getA() const
{
    return A;
}
