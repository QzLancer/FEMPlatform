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

	void addNonlinearMaterial(string name, int bhpoints, double* bdata, double* hdata);
	void addLinearMaterial(string name, double mu);
	virtual void readMeshFile(string meshfile) = 0;
	virtual void createElement2Material() = 0;	/*���õ�Ԫ��������*/
	virtual void bulidGeometry2Load() = 0;		/*���ø���*/
	virtual void buildGeometry2Constrain() = 0;	/*���ñ߽�����*/

protected:
	//stl�е�ָ���������������Ӧ�øĳ�����ָ��
	std::vector<FEMMaterial*> materiallist;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
};

