#include "FEMMeshManager.h"


FEMMeshManager::FEMMeshManager():
	m_num_nodes(0),
	m_num_vtxele(0),
	m_num_edgele(0),
	m_num_triele(0),
	mp_node(nullptr),
	mp_vtxele(nullptr),
	mp_edgele(nullptr),
	mp_triele(nullptr)
{
}

FEMMeshManager::~FEMMeshManager()
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
	//if (model != nullptr) {
	//	gmsh::finalize();
	//}
	gmsh::finalize();
}

void FEMMeshManager::readGeoFile(string geofile)
{
	//获取几何文件拓展名
	string extension;
	for (int i = geofile.size() - 1; i >= 0; --i) {
		if (geofile[i] == '.') {
			extension = geofile.substr(i, geofile.size() - i);
			break;
		}
	}

	if (extension == ".geo") {
		int myargn = 4;
		char* myargv[] = { (char*)"gmsh",(char*)"-format",(char*)"msh2",(char*)"-v",(char*)"0" };
		gmsh::initialize(myargn, myargv);
		gmsh::option::setNumber("General.Terminal", 1);
		gmsh::open(geofile);
 		cout << "Opening model " << geofile << endl << endl;
		printf("Meshing...\n");
		/** 初次分网，重分网必须先有一个网格 **/
		gmsh::model::mesh::generate();
		meshfile = geofile + "_0.msh";
		gmsh::write(meshfile);
	}
}

void FEMMeshManager::meshUnitConvert(double unitratio)
{
	for (int i = 0; i < m_num_nodes; ++i) {
		mp_node[i].x *= unitratio;
		mp_node[i].y *= unitratio;
	}
}

void FEMMeshManager::remesh(double dx, double dy)
{

}

void FEMMeshManager::deleteMeshDomain(int dimension, int domain)
{
	cout << "Delete mesh Domain" << domain <<"...\n";
	pair<int, int> domainpair(dimension, 13);
	gmsh::vectorpair vecpair = { domainpair };
	//gmsh::model::mesh::clear();
}

int FEMMeshManager::getNumofNodes() const
{
	return m_num_nodes;
}

CNode* FEMMeshManager::getNodes() const
{
	return mp_node;
}

int FEMMeshManager::getNumofVtxEle() const
{
	return m_num_vtxele;
}

CVtxElement* FEMMeshManager::getVtxElements() const
{
	return mp_vtxele;
}

int FEMMeshManager::getNumofEdgEle() const
{
	return m_num_edgele;
}

CEdgElement* FEMMeshManager::getEdgElements() const
{
	return mp_edgele;
}

int FEMMeshManager::getNumofTriEle() const
{
	return m_num_triele;
}

CTriElement* FEMMeshManager::getTriElements() const
{
	return mp_triele;
}

int FEMMeshManager::next_int(char** start)
{
	int i;
	char* end;

	i = strtol(*start, &end, 10);
	*start = end;
	return(i);
}
