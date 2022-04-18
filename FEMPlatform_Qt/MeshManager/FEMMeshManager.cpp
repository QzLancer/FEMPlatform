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
	//��ȡ�����ļ���չ��
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
		model = GModel::current();
		/** ���η������ط�����������һ������ **/
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

void FEMMeshManager::remesh(string filename, int current_step, double dx, double dy)
{
	meshfile = "../model/" + filename + "_" + to_string(current_step) + ".msh";
	/** δ����λ�ƾͲ�Ҫ���� **/
	if (fabs(dx) < 1e-10 && fabs(dy) < 1e-10) {
		printf("Skipping remeshing.\n");
		gmsh::write(meshfile);
		return;
	}

	/** ɾ�������ķ��� **/
	GFace* f_xiantie = nullptr;
	GFace* f_air = nullptr;
	for (GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it) {
		if ((*it)->tag() == tag_remesh) {
			f_air = (*it);
		}
		if ((*it)->tag() == tag_armature) {
			f_xiantie = (*it);
		}
	}

	/** �ƶ������ķ��� **/
	if (f_xiantie) {
		printf("moving armature region %d mesh dx=%lf,dy=%lf...\n", tag_armature, dx, dy);
		moveFace(f_xiantie, dx, dy, 0);
	}

	if (f_air) {
		/** ɾ�������ķ��� **/
		printf("deleting air surface %d\n", tag_remesh);
		deleteFaceMesh(f_air);
		/** �Կ��������ط��� **/
		printf("remesh air domain...\n");
		f_air->mesh(true);
	}
	gmsh::write(meshfile);
	readMeshFile(meshfile);
	printf("Finish remesh 2d\n");

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

void FEMMeshManager::deleteFaceMesh(GFace* f)
{
	deMeshGFace dem;
	dem(f);
}

void FEMMeshManager::moveFace(GFace* f, double dx, double dy, double dz) {
	/** �޸��ڲ������ڵ����� **/
	for (std::size_t j = 0; j < f->getNumMeshVertices(); ++j) {
		f->getMeshVertex(j)->setXYZ(f->getMeshVertex(j)->x() + dx,
			f->getMeshVertex(j)->y() + dy,
			f->getMeshVertex(j)->z() + dz);
	}
	/** �޸ı߽����Ϸ����ڵ����� **/
	std::set<MVertex*, MVertexLessThanNum> all_vertices;
	for (auto e : f->edges()) {
		for (auto line : e->lines) {
			MVertex* v1 = line->getVertex(0);
			MVertex* v2 = line->getVertex(1);

			all_vertices.insert(v1);
			all_vertices.insert(v2);
		}
	}
	/** all_vertices��e->getMeshVertex���˶��� **/
	for (std::set<MVertex*, MVertexLessThanNum>::iterator ite = all_vertices.begin(); ite != all_vertices.end(); ite++) {
		(*ite)->setXYZ((*ite)->x() + dx, (*ite)->y() + dy, (*ite)->z() + dz);
	}
	/** �޸ļ��ζ������� **/
	for (auto v : f->vertices()) {
		for (std::size_t j = 0; j < v->getNumMeshVertices(); ++j) {
			GPoint p(v->x() + dx, v->y() + dy, v->z() + dz);
			v->setPosition(p);
		}
	}
	/** ����distance field **/
	updateField();
}

void FEMMeshManager::updateField() {
	/** ����distance field **/
	FieldManager* fM = GModel::current()->getFields();
	std::map<int, Field*>::iterator iter;
	for (iter = fM->begin(); iter != fM->end(); ++iter) {
		iter->second->update_needed = true;
		iter->second->update();
	}
}

void deleteFaceMesh(GFace* f) {
	deMeshGFace dem;
	dem(f);
}