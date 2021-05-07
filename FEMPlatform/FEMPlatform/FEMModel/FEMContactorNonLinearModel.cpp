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
	meshfile = "../../model/model3656.mphtxt";
}
 
void FEMContactorNonLinearModel::createElement2Material()
{
    //添加材料，后续应该是在用户界面添加，然后一直以文本的形式保存
    addLinearMaterial("Air", 4 * PI * 1e-7);
    double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4 , 3};
    double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3,
        20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 , 1000000};
    addNonlinearMaterial("Soft Iron1", 16, bdata, hdata);

    //double* hdata = new double[] {  0, 10, 100, 150, 200, 250, 300, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 900, 1900, 2900, 3900, 4900, 5900, 6900, 7900, 8900, 9900  };
    //double* bdata = new double[] { 0, 0.0599718475592522, 0.597090304698867, 0.891647197655098, 1.17932987869717, 1.44849641137579, 1.66011089327454, 1.76457710473790, 1.77545758128628, 1.78436476327606, 1.79172378239890, 1.79786552945965, 1.80304417672406, 1.80745433989896, 1.81124544595357, 1.81453294220234, 1.81740674973144, 1.81993754230050, 1.82218138737627, 1.82418318028023, 1.82597919696382, 1.82759900417761, 1.82906689986371, 1.83040300831993, 1.83162411995517, 1.83274434065704, 1.83377559810748, 1.83472803974086, 1.83561034796116, 1.83642999167750, 1.83719342844891, 1.83790626803611, 1.83857340558121, 1.83919913072030, 1.83978721749903, 1.84034099887956, 1.84086342880472, 1.84135713415656, 1.84182445846260, 1.84226749882734, 1.84268813727457, 1.84308806745612, 1.84346881750189, 1.84383176964217, 1.84417817711895, 1.84450917881082, 1.84482581192233, 1.84512902302871, 1.84541967771808, 1.84569856903387, 1.84596642488713, 1.84622391458191, 1.84647165457468, 1.84851935249073, 1.85578788972131, 1.85843034336099, 1.86030242999031, 1.86190677974238, 1.86338688482529, 1.86479924261073, 1.86617060520324, 1.86751527801887, 1.86884160416511  };
    //addNonlinearMaterial("Soft Iron", 63, bdata, hdata);


    //针对当前模型设置材料，和上面的添加材料功能应该分离到两个不同的类
    materialmap[1] = materiallist[0];   //空气
    materialmap[2] = materiallist[0];   //压缩空气
    materialmap[3] = materiallist[1];   //衔铁
    materialmap[4] = materiallist[1];   //轭铁
    materialmap[5] = materiallist[0];   //线圈
}

void FEMContactorNonLinearModel::bulidGeometry2Load()
{
    loadmap[5] = 8e6;
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
    FEMMovingPart moving1;
    moving1.direction[1] = 1;
    movingmap[3] = moving1;
}
