#pragma once
#include "FEMDataType.h"
#include "FEMMaterial.h"
#include "FEMBoundary.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

class FEMGenerator
{
public:
	FEMGenerator();
	~FEMGenerator();
	void addNonlinearMaterial(string name, int bhpoints, double* bdata, double* hdata);
	void addLinearMaterial(string name, double mu);
	virtual void readMeshFile(string meshfile) = 0;
	virtual void createElement2Material() = 0;	/*���õ�Ԫ��������*/
	virtual void bulidGeometry2Load() = 0;		/*���ø���*/
	virtual void buildGeometry2Constrain() = 0;	/*���ñ߽�����*/

	int getNumofNodes() const;
	CNode* getNodes() const;
	int getNumofVtxEle() const;
	CVtxElement* getVtxElements() const;
	int getNumofEdgEle() const;
	CEdgElement* getEdgElements() const;
	int getNumofTriEle() const;
	CTriElement* getTriElements() const;
	std::map<int, FEMMaterial*> getMaterial() const;
	std::map<int, double> getLoad() const;
	std::map<int, FEMBoundary*> getBoundary() const;

protected:
	//ָ��������յĻ����Ǻܺã�����Ӧ�øĳ�����ָ���STL
	//��Ԫ������Ҫ��ģ�������Ԫ����������Builderģʽ����Ԫ����
	int m_num_nodes;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;
	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;
	std::vector<FEMMaterial*> materiallist;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
};

