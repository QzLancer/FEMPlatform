#include "FEMContactorLinearModel.h"

void FEMContactorLinearModel::setModelName()
{
    modelname = "ContactLinear";
}

void FEMContactorLinearModel::setDimension()
{
    dimension = DIMENSION::TWO;
}

void FEMContactorLinearModel::setFile()
{
    meshfile = "../../model/model.mphtxt";
}

void FEMContactorLinearModel::createElement2Material()
{
    //��Ӳ��ϣ�����Ӧ�������û�������ӣ�Ȼ��һֱ���ı�����ʽ����
    addLinearMaterial("Air", 4 * PI * 1e-7);
    //double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4 };
    //double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3,
    //    20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
    //addNonlinearMaterial("Soft Iron", 16, bdata, hdata);
    addLinearMaterial("Soft Iron", 1. / 664.);  //���ͳ���/���ͳ����õ��Ļ�������

    //��Ե�ǰģ�����ò��ϣ����������Ӳ��Ϲ���Ӧ�÷��뵽������ͬ����
    materialmap[1] = materiallist[0];
    materialmap[2] = materiallist[0];
    materialmap[3] = materiallist[1];
    materialmap[4] = materiallist[1];
    materialmap[5] = materiallist[0];
}

void FEMContactorLinearModel::bulidGeometry2Load()
{
    loadmap[5] = 8e5;
}

void FEMContactorLinearModel::buildGeometry2Constrain()
{
    boundarymap[1] = new FEMAxialSymmetry;
    boundarymap[3] = new FEMAxialSymmetry;
    boundarymap[5] = new FEMAxialSymmetry;
    std::vector<int> id = { 2, 31, 32, 37, 38, 48, 49, 52, 53, 59, 60, 61, 62, 63, 64 };   //��ô��֤�߽�����ȷ�ģ�
    for (auto a : id) {
        boundarymap[a] = new FEMMagneticInsulation;
    }
}

void FEMContactorLinearModel::setUnitRatio()
{
    unitratio = 1;
}
