//#include <QCoreApplication>
#include "FEMCore.h"
#include "FEMModel/FEMContactorLinearModel.h"
#include "FEMModel/FEMContactorNonLinearModel.h"
#include "FEMModel/FEMRelay1250Model.h"
#include "FEMModel/FEMRelay1250LinearModel.h"
#include "FEMModel/RelayDynamicModel.h"
#include "FEMModel/FEMTrans3PhaseModel.h"
#include "FEMModel/RelayModelwithoutBand.h"

#include <iostream>
#include <time.h>
#include <direct.h>

std::string analysistype = "static";
std::string solvestrategy = "DD";
std::string matrixsolver = "SuperLU_MT";

int maxitersteps = 10000;
double maxerror = 1e-6;

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    /*
    Model:
        FEMContactorLinearModel
        FEMContactorNonLinearModel
        RelayModelwithoutBand
        FEMRelay1250Model
        RelayDynamicModel
        RelayModelwithoutBand
    */
    // SoftIron1 对应的是师兄的BH曲线

    char* buffer = getcwd(NULL, 0);
    printf("%s\n", buffer);
    free(buffer);

    FEMCore core;
    FEMModel* model = new FEMContactorLinearModel;
    model->init();

    core.setModel(model);
    core.setAnalysisType(analysistype);
    core.setFEMSolveStrategy(solvestrategy);
    core.setMaxIterSteps(maxitersteps);
    core.setMatrixSolver(matrixsolver);
    core.setMaxIterSteps(maxitersteps);
    core.setMaxError(maxerror);

    core.solve(analysistype);


    core.postprocess();

    return 0;
}
