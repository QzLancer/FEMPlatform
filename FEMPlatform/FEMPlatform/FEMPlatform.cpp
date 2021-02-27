// FemPlatform.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "FEMCore.h"

std::string dimension = "2D";
std::string meshfile = "../../model/model.mphtxt";  //文件的起始目录是.sln文件？？
std::string e2mfile;    //可以考虑通过文件的形式添加材料、负载、约束
std::string g2lfile;
std::string g2cfile;
std::string analysistype = "static";
std::string solvestrategy = "NR";
std::string matrixsolver = "SuperLU_MT";

int maxitersteps = 200;

int main()
{
    FEMCore core;
    core.setDimension(dimension);
    core.readMeshData(meshfile);
    core.createElement2Material();
    core.bulidGeometry2Load();
    core.buildGeometry2Constrain();
    core.setAnalysisType(analysistype);
    core.setFEMSolveStrategy(solvestrategy);
    core.setMaxIterSteps(maxitersteps);
    core.setMatrixSolver(matrixsolver);
    core.solve();
    core.postoperation();
}
