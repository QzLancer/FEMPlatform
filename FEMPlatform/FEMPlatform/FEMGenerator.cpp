#include "FEMGenerator.h"


FEMGenerator::FEMGenerator():
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

FEMGenerator::~FEMGenerator()
{
	if (mp_triele != nullptr)
		delete[] mp_triele;
	if (mp_edgele != nullptr)
		delete[] mp_edgele;
	if (mp_vtxele != nullptr)
		delete[] mp_vtxele;
	if (mp_node != nullptr)
		delete[] mp_node;
}

void FEMGenerator::addNonlinearMaterial(string name, int bhpoints, double* bdata, double* hdata)
{
	FEMMaterial* material = new FEMMaterial;
	material->setName(name);
	material->setLinearFlag(false);
	material->setBHpoints(bhpoints);
	material->setBdata(bdata);
	material->setHdata(hdata);
	materiallist.push_back(material);
}

void FEMGenerator::addLinearMaterial(string name, double mu)
{
	FEMMaterial* material = new FEMMaterial;
	material->setBHpoints(1);
	material->setLinearFlag(true);
	material->setmu(mu);
	materiallist.push_back(material);
}

int FEMGenerator::getNumofNodes() const
{
	return m_num_nodes;
}

CNode* FEMGenerator::getNodes() const
{
	return mp_node;
}

int FEMGenerator::getNumofVtxEle() const
{
	return m_num_vtxele;
}

CVtxElement* FEMGenerator::getVtxElements() const
{
	return mp_vtxele;
}

int FEMGenerator::getNumofEdgEle() const
{
	return m_num_edgele;
}

CEdgElement* FEMGenerator::getEdgElements() const
{
	return mp_edgele;
}

int FEMGenerator::getNumofTriEle() const
{
	return m_num_triele;
}

CTriElement* FEMGenerator::getTriElements() const
{
	return mp_triele;
}

std::map<int, FEMMaterial*> FEMGenerator::getMaterial() const
{
	return materialmap;
}

std::map<int, double> FEMGenerator::getLoad() const
{
	return loadmap;
}

std::map<int, FEMBoundary*> FEMGenerator::getBoundary() const
{
	return boundarymap;
}
