#pragma once
#include "../FEMDataType.h"
#include "../FEMMaterial.h"
#include "../FEMBoundary.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

class FEMMeshManager
{
public:
	FEMMeshManager();
	~FEMMeshManager();
	virtual void readMeshFile(string meshfile) = 0;


	int getNumofNodes() const;
	CNode* getNodes() const;
	int getNumofVtxEle() const;
	CVtxElement* getVtxElements() const;
	int getNumofEdgEle() const;
	CEdgElement* getEdgElements() const;
	int getNumofTriEle() const;
	CTriElement* getTriElements() const;

protected:
	//指针这块掌握的还不是很好，后续应该改成智能指针和STL
	//单元部分需要大改，创建单元抽象类，Builder模式处理单元类型
	int m_num_nodes;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;
	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;
};

