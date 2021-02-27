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
	//ָ��������յĻ����Ǻܺã�����Ӧ�øĳ�����ָ���STL
	//��Ԫ������Ҫ��ģ�������Ԫ�����࣬Builderģʽ����Ԫ����
	int m_num_nodes;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;
	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;
};

