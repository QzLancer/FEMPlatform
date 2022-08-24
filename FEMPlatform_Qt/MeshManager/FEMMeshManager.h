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

	/// @brief 读取并解析.geo格式的几何文件，通过Gmsh初次分网，生成.msh网格文件
	/// @param geofile .geo格式几何文件路径
	virtual void readGeoFile(string geofile);

	/// @brief 读取.mphtxt格式或.msh的网格文件
	/// @param meshfile 网格文件路径
	virtual void readMeshFile(string meshfile = "") = 0;

	/// @brief 网格单位转换
	/// @param unitratio 转换系数
	virtual void meshUnitConvert(double unitratio);

	/// @brief 重分网，该函数目前还不完善，需要进一步设计（例如输入参数应该包含需要重分网的区域的编号和位移区域的编号）
	/// @param filename 重分网后保存的网格文件名称
	/// @param current_step 当前时间部
	/// @param dx x方向位移
	/// @param dy y方向位移
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
	//指针这块掌握的还不是很好，后续应该改成智能指针和STL
	//单元部分需要大改，创建单元抽象类，Builder模式处理单元类型
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

	/// @brief 指定remesh区域编号，建议将该变量修改为remesh函数的一个参数
	int tag_remesh{5};
	/// @brief 指定armature运动区域编号，建议将该变量修改为remesh函数的一个参数
	int tag_armature{1};
};

