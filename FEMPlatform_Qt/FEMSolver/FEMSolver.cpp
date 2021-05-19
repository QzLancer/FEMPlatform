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
    maxitersteps(100),
	current(nullptr),
	dis(nullptr),
	velocity(nullptr),
	acc(nullptr),
	magneticforce(nullptr),
	springforce(nullptr),
	mass(0.024)
{

}

FEMSolver::~FEMSolver()
{
	if (mp_triele != nullptr) {
		delete[] mp_triele;
		mp_triele = nullptr;
	}
	if (mp_edgele != nullptr) {
		delete[] mp_edgele;
		mp_edgele = nullptr;
	}
	if (mp_vtxele != nullptr) {
		delete[] mp_vtxele;
		mp_vtxele = nullptr;
	}
	if (mp_node != nullptr) {
		delete[] mp_node;
		mp_node = nullptr;
	}
	if (strategy != nullptr) {
		delete strategy;
		strategy = nullptr;
	}
	if (matsolver != nullptr) {
		delete matsolver;
		matsolver = nullptr;
	}
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
	cout << "FEMSolver::setNodes\n";
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
	cout << "FEMSolver::setTriElements\n";
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

void FEMSolver::setMaxError(const double _error)
{
	maxerror = _error;
}

void FEMSolver::writeVtkFile(std::string _name)
{
	writeGeometryVtkFile(_name);
	std::string name = std::string("D:/femplatform/result/");
	name += _name;
	name += "_";
	name += /*to_string(m_num_nodes) +*/ std::string(".vtk");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: opening file!" << endl;
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
		fprintf(fp, "%f\n", mp_node[i].A);
	}

	fprintf(fp, "SCALARS F double 1\n");
	fprintf(fp, "LOOKUP_TABLE %s\n", "Ftable");
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%f\n", mp_node[i].NodeForce);
	}

	fprintf(fp, "\nVECTORS %s double\n", "Fvector");
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%lf %lf %lf\n", mp_node[i].NodeForcex, mp_node[i].NodeForcey, 0.0);
	}

	//单元标量磁感应强度
	fprintf(fp, "\nCELL_DATA %d\n", m_num_triele);
	fprintf(fp, "SCALARS %s double %d\n", "Bnorm", 1);
	fprintf(fp, "LOOKUP_TABLE %s\n", "Btable");
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		fprintf(fp, "%lf\n", mp_triele[i_tri].B);
	}

	//单元矢量磁感应强度
	fprintf(fp, "\nVECTORS %s double\n", "Bvector");
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		fprintf(fp, "%lf %lf %lf\n", mp_triele[i_tri].Bx, mp_triele[i_tri].By, 0.0);
	}

	//单元磁导率
	fprintf(fp, "SCALARS %s double %d\n", "Mu", 1);
	fprintf(fp, "LOOKUP_TABLE %s\n", "Mutable");
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		fprintf(fp, "%.10f\n", mp_triele[i_tri].material->getMu(mp_triele[i_tri].B));
	}
	fclose(fp);
}

void FEMSolver::writeTxtFile(std::string _name)
{
	std::string name = std::string("D:/femplatform/result/");
	name += _name;
	name += "_";
	name += to_string(m_num_nodes) + std::string(".txt");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: opening file!" << endl;
		exit(0);
	}
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%lf %lf %lf\n", mp_node[i].x, mp_node[i].y, mp_node[i].A);
	}
}

void FEMSolver::writeVtkFileNoAir(std::string _name, vector<int> air_domain)
{
	writeGeometryVtkFile(_name);
	std::string name = std::string("D:/femplatform/result/");
	name += _name;
	name += "_noAir";
	name += /*to_string(m_num_nodes) +*/ std::string(".vtk");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: opening file!" << endl;
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

	/** 计算非空气单元数目 **/
	int num_cell;
	bool airflag = false;
	vector<int> eleid;
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		airflag = false;
		for (auto a : air_domain) {
			if (mp_triele[i_tri].domain == a) {
				airflag = true;
				break;
			}
		}
		if (airflag == false) {
			eleid.push_back(i_tri);
		}
	}
	num_cell = eleid.size();

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
	fprintf(fp, "\nCELLS %d %d\n", num_cell, 4 * num_cell);
	for (int i = 0; i < num_cell; ++i) {
		fprintf(fp, "3 %d %d %d\n", mp_triele[eleid[i]].n[0], mp_triele[eleid[i]].n[1], mp_triele[eleid[i]].n[2]);
	}
	fprintf(fp, "\nCELL_TYPES %d\n", num_cell);
	int type = 5;
	for (int i = 0; i < num_cell; ++i) {
		fprintf(fp, "%d\n", type);
	}
	//节点磁矢位
	fprintf(fp, "\nPOINT_DATA %d\n", m_num_nodes);
	fprintf(fp, "SCALARS A double 1\n");
	fprintf(fp, "LOOKUP_TABLE %s\n", "Atable");
	for (int i = 0; i < m_num_nodes; ++i) {
		fprintf(fp, "%f\n", mp_node[i].A);
	}

	//单元标量磁感应强度
	fprintf(fp, "\nCELL_DATA %d\n", num_cell);
	fprintf(fp, "SCALARS %s double %d\n", "Bnorm", 1);
	fprintf(fp, "LOOKUP_TABLE %s\n", "Btable");
	for (int i_tri = 0; i_tri < num_cell; ++i_tri) {
		fprintf(fp, "%lf\n", mp_triele[eleid[i_tri]].B);
	}

	fprintf(fp, "\nVECTORS %s double\n", "Bvector");
	for (int i_tri = 0; i_tri < num_cell; ++i_tri) {
		fprintf(fp, "%lf %lf %lf\n", mp_triele[eleid[i_tri]].Bx, mp_triele[eleid[i_tri]].By, 0.0);
	}
	fclose(fp);
}

void FEMSolver::writeGeometryVtkFile(std::string _name)
{
	std::string name = std::string("D:/femplatform/result/");
	name += _name;
	name += "_";
	//name += to_string(m_num_nodes);
	name += "Geo" + std::string(".vtk");
	FILE* fp = nullptr;
	int err;
	char ch[256];
	err = fopen_s(&fp, name.c_str(), "w");
	if (!fp) {
		std::cout << "Error: opening file!" << endl;
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
	//线单元
	fprintf(fp, "\nCELLS %d %d\n", m_num_edgele, 3 * m_num_edgele);
	for (int i = 0; i < m_num_edgele; ++i) {
		fprintf(fp, "2 %d %d\n", mp_edgele[i].n[0], mp_edgele[i].n[1]);
	}
	fprintf(fp, "\nCELL_TYPES %d\n", m_num_edgele);
	int type = 3;
	for (int i = 0; i < m_num_edgele; ++i) {
		fprintf(fp, "%d\n", type);
	}
	fclose(fp);
}

void FEMSolver::setDeformedDomain(const std::vector<int> _deformedlist)
{
	deformedlist = _deformedlist;
}

void FEMSolver::setMovingPart(const std::map<int, FEMMovingPart*> _movingmap)
{
	movingmap = _movingmap;
}

void FEMSolver::setMeshManager(FEMMeshManager* _meshmanager)
{
	meshmanager = _meshmanager;
}

void FEMSolver::updateLoadmap(int domain, double current)
{
	double Jor = current * materialmap[domain]->getFEMCoil().Nc / materialmap[domain]->getFEMCoil().width / materialmap[domain]->getFEMCoil().height;
	loadmap[domain] = Jor;
}

//std::vector<double> FEMSolver::getA() const
//{
//    return A;
//}
