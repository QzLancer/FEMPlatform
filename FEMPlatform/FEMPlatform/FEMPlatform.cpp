﻿// FemPlatform.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "FEMCore.h"
#include "FEMModel/FEMContactorLinearModel.h"
#include "FEMModel/FEMContactorNonLinearModel.h"
#include "FEMModel/FEMRelay1250Model.h"
#include "FEMModel/FEMRelay1250LinearModel.h"

#include <iostream>
#include <time.h>

std::string analysistype = "static";
std::string solvestrategy = "NDDR";
std::string matrixsolver = "SuperLU_MT";

int maxitersteps = 20000;
double maxerror = 1e-7;

int main()
{
    FEMCore core;
    FEMModel* model = new FEMRelay1250Model;
    model->init();

    core.setModel(model);
    core.setAnalysisType(analysistype);
    core.setFEMSolveStrategy(solvestrategy);
    core.setMaxIterSteps(maxitersteps);
    core.setMatrixSolver(matrixsolver);
    core.setMaxIterSteps(maxitersteps);
    core.setMaxError(maxerror);
    clock_t start, end;
    start = clock();
    core.solve();
    end = clock();
    cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
    core.postprocess();
}
