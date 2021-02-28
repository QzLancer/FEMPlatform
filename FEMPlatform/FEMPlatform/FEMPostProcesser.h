#pragma once

#include "MeshManager/FEMMeshManager.h"
#include "FEMDataType.h"

#include <string>
#include <vector>
#include <iostream>

class FEMPostProcesser
{
public:
	//��ά��д���ʸλA
	void writeVtkFile(std::string _name, FEMMeshManager* _meshmanager, std::vector<double> _A);
};

