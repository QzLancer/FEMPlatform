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

	/// @brief ��ȡ������.geo��ʽ�ļ����ļ���ͨ��Gmsh���η���������.msh�����ļ�
	/// @param geofile .geo��ʽ�����ļ�·��
	virtual void readGeoFile(string geofile);

	/// @brief ��ȡ.mphtxt��ʽ��.msh�������ļ�
	/// @param meshfile �����ļ�·��
	virtual void readMeshFile(string meshfile = "") = 0;

	/// @brief ����λת��
	/// @param unitratio ת��ϵ��
	virtual void meshUnitConvert(double unitratio);

	/// @brief �ط������ú���Ŀǰ�������ƣ���Ҫ��һ����ƣ������������Ӧ�ð�����Ҫ�ط���������ı�ź�λ������ı�ţ�
	/// @param filename �ط����󱣴�������ļ�����
	/// @param current_step ��ǰʱ�䲿
	/// @param dx x����λ��
	/// @param dy y����λ��
	virtual void remesh(string filename, int current_step, double dx, double dy);


	int getNumofNodes() const;
	CNode* getNodes() const;
	int getNumofVtxEle() const;
	CVtxElement* getVtxElements() const;
	int getNumofEdgEle() const;
	CEdgElement* getEdgElements() const;
	int getNumofTriEle() const;
	CTriElement* getTriElements() const;


private:	
	void deleteFaceMesh(GFace* f);
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
	GModel* model;

	/// @brief ָ��remesh�����ţ����齫�ñ����޸�Ϊremesh������һ������
	int tag_remesh{5};
	/// @brief ָ��armature�˶������ţ����齫�ñ����޸�Ϊremesh������һ������
	int tag_armature{1};
};

