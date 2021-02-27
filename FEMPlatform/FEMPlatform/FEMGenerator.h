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
	virtual void createElement2Material() = 0;	/*设置单元材料类型*/
	virtual void bulidGeometry2Load() = 0;		/*设置负载*/
	virtual void buildGeometry2Constrain() = 0;	/*设置边界条件*/

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
	//指针这块掌握的还不是很好，后续应该改成智能指针和STL
	//单元部分需要大改，创建单元抽象类来，Builder模式处理单元类型
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

