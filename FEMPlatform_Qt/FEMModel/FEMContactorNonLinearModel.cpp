#include "FEMContactorNonLinearModel.h"

void FEMContactorNonLinearModel::setModelName()
{
	modelname = "ContactnonLinear";
}

void FEMContactorNonLinearModel::setDimension()
{
	dimension = DIMENSION::D2AXISM;
}

void FEMContactorNonLinearModel::setFile()
{
	meshfile = "../model/model8769.mphtxt";
}
 
void FEMContactorNonLinearModel::createElement2Material()
{
    //添加材料，后续应该是在用户界面添加，然后一直以文本的形式保存
    addLinearMaterial("Air", 4 * PI * 1e-7);
    double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
    double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
    addNonlinearMaterial("Soft Iron", 16, bdata, hdata);

    //针对当前模型设置材料，和上面的添加材料功能应该分离到两个不同的类
    materialmap[1] = materiallist[0];   //空气
    materialmap[2] = materiallist[0];   //压缩空气
    materialmap[3] = materiallist[1];   //衔铁
    materialmap[4] = materiallist[1];   //轭铁
    materialmap[5] = materiallist[0];   //线圈
}

void FEMContactorNonLinearModel::bulidGeometry2Load()
{
    loadmap[5] = 8372903.023;
    //loadmap[5] = 14372903.023;
}

void FEMContactorNonLinearModel::buildGeometry2Constrain()
{
    boundarymap[1] = new FEMAxialSymmetry;
    boundarymap[3] = new FEMAxialSymmetry;
    boundarymap[5] = new FEMAxialSymmetry;
    std::vector<int> id = { 2, 7, 31, 32, 37, 38, 48, 49, 52, 53, 59, 60, 61, 62, 63, 64 };   //怎么验证边界是正确的？
    for (auto a : id) {
        boundarymap[a] = new FEMMagneticInsulation;
    }
}

void FEMContactorNonLinearModel::setUnitRatio()
{
    unitratio = 1;
}

void FEMContactorNonLinearModel::buildGeometry2Deformed()
{
    deformedlist.push_back(2);
}

void FEMContactorNonLinearModel::buildGeometry2MovingPart()
{
    FEMMovingPart* moving1 = new FEMMovingPart;
    moving1->direction[1] = 1;
    movingmap[3] = moving1;
}
