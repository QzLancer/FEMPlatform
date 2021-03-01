// FemPlatform.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "FEMCore.h"
#include "FEMModel/FEMContactorLinearModel.h"
#include "FEMModel/FEMContactorNonLinearModel.h"

std::string analysistype = "static";
std::string solvestrategy = "NR";
std::string matrixsolver = "SuperLU_MT";

int maxitersteps = 200;

int main()
{
    FEMCore core;
    FEMModel* model = new FEMContactorLinearModel;
    model->init();

    core.setModel(model);
    core.setAnalysisType(analysistype);
    core.setFEMSolveStrategy(solvestrategy);
    core.setMaxIterSteps(maxitersteps);
    core.setMatrixSolver(matrixsolver);
    core.setMaxIterSteps(200);
    core.solve();
    core.postprocess();
}
