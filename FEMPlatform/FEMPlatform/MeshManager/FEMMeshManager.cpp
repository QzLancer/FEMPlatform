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
	if (mp_triele != nullptr)
		delete[] mp_triele;
	if (mp_edgele != nullptr)
		delete[] mp_edgele;
	if (mp_vtxele != nullptr)
		delete[] mp_vtxele;
	if (mp_node != nullptr)
		delete[] mp_node;
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