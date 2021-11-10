#include <QCoreApplication>
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

std::string analysistype = "static";
std::string solvestrategy = "NDDR";
std::string matrixsolver = "SuperLU_MT";

int maxitersteps = 100000;
double maxerror = 1e-5;

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    /*
    Model:
        FEMContactorLinearModel
        FEMContactorNonLinearModel
        RelayModelwithoutBand
    */

    FEMCore core;
    FEMModel* model = new FEMContactorNonLinearModel;
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
