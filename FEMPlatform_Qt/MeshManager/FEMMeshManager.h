#pragma once
#include "../FEMDataType.h"
#include "../FEMMaterial.h"
#include "../FEMBoundary.h"
#include "meshGFace.h"
#include "meshGRegionNetgen.h"
#include "meshGRegion.h"
#include "Generator.h"
#include "GModel.h"
#include "gmsh.h"
#include "MLine.h"
#include "Field.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>

using namespace std;

class FEMMeshManager
{
public:
	FEMMeshManager();
	virtual ~FEMMeshManager();
	virtual void readGeoFile(string geofile);
	virtual void readMeshFile(string meshfile = "") = 0;
	virtual void meshUnitConvert(double unitratio);
	virtual void remesh(double dx, double dy);
	void deleteMeshDomain(int dimension, int domain);

	int getNumofNodes() const;
	CNode* getNodes() const;
	int getNumofVtxEle() const;
	CVtxElement* getVtxElements() const;
	int getNumofEdgEle() const;
	CEdgElement* getEdgElements() const;
	int getNumofTriEle() const;
	CTriElement* getTriElements() const;


private:
	void moveFace(GFace* f, double dx, double dy, double dz);
	void updateField();

protected:
	//ָ��������յĻ����Ǻܺã�����Ӧ�øĳ�����ָ���STL
	//��Ԫ������Ҫ��ģ�������Ԫ�����࣬Builderģʽ����Ԫ����
	int m_num_nodes;
	int m_num_ele;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;
	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;

	string meshfile;

	int next_int(char **start);
	static GModel* model;
};

