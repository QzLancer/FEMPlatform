#include "FEMRelay1250LinearModel.h"

void FEMRelay1250LinearModel::setModelName()
{
	modelname = "JRS1250Linear";
}

void FEMRelay1250LinearModel::setDimension()
{
	dimension = DIMENSION::D2AXISM;
}

void FEMRelay1250LinearModel::setFile()
{
	//geofile = "../../model/geo/JRS1250/JRS1250bgmband.geo";
    meshfile = "../../model/geo/JRS1250/JRS1250bgmband.geo_0_22888.msh";
}

void FEMRelay1250LinearModel::createElement2Material()
{
    addLinearMaterial("Air", 4 * PI * 1e-7);
    //double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4 };
    //double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3,
    //    20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
    //addLinearMaterial("PM", 1.12 * 4 * PI * 1e-7, 1.3 / 1.12 / (4 * PI * 1e-7), 0);
    addPermMagentMaterial("PM", 1.12 * (4 * PI * 1e-7), 1.3 / 1.12 / (4 * PI * 1e-7), 1.21, PI);
    addLinearMaterial("Soft Iron", 1. / 664.);  //整型常数/整型常数得到的还是整型
    FEMCoil coil;
    coil.xyz = 1;
    coil.direction = 1;
    coil.Nc = 60;
    coil.Rc = 3.9;
    coil.Vc = 27;
    coil.type = CICULARCOIL;
    coil.height = 17e-3;
    coil.center[0] = 0;
    coil.center[1] = 0;
    coil.center[2] = 0;
    coil.r0 = 0;
    coil.width = 3e-3;
    coil.Jor = coil.Vc / coil.Rc * coil.Nc / coil.width / coil.height;
    coil.tau = coil.Nc / coil.width / coil.height;
    addCoil("Coil", coil);


    materialmap[1] = materiallist[2];   //衔铁
    materialmap[2] = materiallist[0];   //线圈
    materialmap[3] = materiallist[2];   //铁芯
    materialmap[4] = materiallist[2];   //铁芯
    materialmap[5] = materiallist[2];   //铁芯
    materialmap[6] = materiallist[1];   //永磁
    materialmap[7] = materiallist[2];   //铁芯（导磁环）
    materialmap[8] = materiallist[2];   //铁芯
    materialmap[9] = materiallist[0];   //外部空气
    materialmap[10] = materiallist[0];  //空气
    materialmap[11] = materiallist[0];  //空气
    materialmap[12] = materiallist[0];  //空气
    materialmap[13] = materiallist[0];  //可压缩空气
}

void FEMRelay1250LinearModel::bulidGeometry2Load()
{
    //for (auto iter = materialmap.begin(); iter != materialmap.end(); ++iter) {
    //    loadmap[iter->first] = iter->second->getFEMCoil().Jor;
    //}
    loadmap[2] = 8144796.38;
}

void FEMRelay1250LinearModel::buildGeometry2Constrain()
{
    std::vector<int> id = { 48, 52, 53, 61, 62 };
    for (auto a : id) {
        boundarymap[a] = new FEMAxialSymmetry;
    }
    id = { 49, 50 };
    for (auto a : id) {
        boundarymap[a] = new FEMMagneticInsulation;
    }
}

void FEMRelay1250LinearModel::setUnitRatio()
{
    unitratio = 0.001;
}

void FEMRelay1250LinearModel::buildGeometry2Deformed()
{
    deformedlist.push_back(13);
}

void FEMRelay1250LinearModel::buildGeometry2MovingPart()
{
    FEMMovingPart* moving1 = new FEMMovingPart;
    moving1->direction[1] = 1;
    moving1->limit[1].min = -0.0062;
    movingmap[1] = moving1;
}
