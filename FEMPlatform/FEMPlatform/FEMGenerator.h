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
	virtual void createElement2Material() = 0;	/*设置单元材料类型*/
	virtual void bulidGeometry2Load() = 0;		/*设置负载*/
	virtual void buildGeometry2Constrain() = 0;	/*设置边界条件*/

protected:
	//stl中的指针如何析构？后续应该改成智能指针
	std::vector<FEMMaterial*> materiallist;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
};

