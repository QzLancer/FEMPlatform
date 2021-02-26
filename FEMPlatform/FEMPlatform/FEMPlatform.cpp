// NodalDomainDecomp.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "FEMCore.h"

std::string dimension = "2D";
std::string meshfile = "../../model/model.mphtxt";
std::string e2mfile;
std::string g2lfile;
std::string g2cfile;
std::string analysistype = "ElectromagneticStatic";
std::string solver = "NR";

int main()
{
    FEMCore core;
    core.setDimension(dimension);
    core.readMeshData(meshfile);
    core.createElement2Material();
    core.bulidGeometry2Load();
    core.buildGeometry2Constrain();
    core.setAnalysisType(analysistype);
    core.setFEMSolver(solver);
    core.postoperation();
}
