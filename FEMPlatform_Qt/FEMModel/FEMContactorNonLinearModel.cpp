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
	meshfile = "../model/model14120.mphtxt";
}
 
void FEMContactorNonLinearModel::createElement2Material()
{
    //添加材料，后续应该是在用户界面添加，然后一直以文本的形式保存
    addLinearMaterial("Air", 4 * PI * 1e-7);
    double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
    double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
    //double* bdata = new double[] {0.000000, 0.222067, 0.434933, 0.629400, 0.796267, 0.929125, 1.032733, 1.114642, 1.182400, 1.242315, 1.295717, 1.342694, 1.383333, 1.417839, 1.446875, 1.471224, 1.491667, 1.508984, 1.523958, 1.537370, 1.550000, 1.562500, 1.575000, 1.587500, 1.600000, 1.612500, 1.625000, 1.637500, 1.650000, 1.662500, 1.675000, 1.687500, 1.700000, 1.712500, 1.725000, 1.737500, 1.750000, 1.762500, 1.775000, 1.787500, 1.800000, 1.812500, 1.825000, 1.837500, 1.850000, 1.862500, 1.875000, 1.887500, 1.900000, 1.912500, 1.925000, 1.937500, 1.950000, 1.962500, 1.975000, 1.987500, 2.000000, 2.012500, 2.025000, 2.037500, 2.050000, 2.062435, 2.074479, 2.085742, 2.095833, 2.104492, 2.111979, 2.118685, 2.125000, 2.131250, 2.137500, 2.143750, 2.150000, 2.156250, 2.162500, 2.168750, 2.175000, 2.181315, 2.188021, 2.195508, 2.204167, 2.214205, 2.225098, 2.236138, 2.246617, 2.255991, 2.264375, 2.272046, 2.279283, 2.286351, 2.293460, 2.300811, 2.308600, 2.317014, 2.326185, 2.336235, 2.347283, 2.359476, 2.373065, 2.388325, 2.405533, 2.424779, 2.445402, 2.466553, 2.487383, 2.507160, 2.525606, 2.542560, 2.557862, 2.573476, 2.599873, 2.649647, 2.735397, 2.865594, 3.032219, 3.223130, 3.426183};
    //double* hdata = new double[] {0.000000, 47.022105, 92.853387, 136.303027, 176.180200, 211.831047, 244.749548, 276.966646, 310.513283, 347.238179, 388.261158, 434.519826, 486.951783, 546.583123, 614.793885, 693.052601, 782.827800, 885.230686, 999.943162, 1126.289808, 1263.595200, 1411.470874, 1570.676190, 1742.257465, 1927.261017, 2127.462688, 2347.556421, 2592.965685, 2869.113950, 3180.557626, 3528.384892, 3912.816870, 4334.074683, 4790.601011, 5273.724767, 5772.996422, 6277.966450, 6780.824385, 7284.316015, 7793.826191, 8314.739767, 8852.425094, 9412.184533, 9999.303948, 10619.069200, 11276.632610, 11976.612340, 12723.493020, 13521.759250, 14374.508630, 15279.290610, 16232.267620, 17229.602070, 18268.616200, 19351.271590, 20480.689600, 21659.991630, 22891.216730, 24172.074580, 25499.192550, 26869.198000, 28276.846980, 29709.410340, 31152.287590, 32590.878280, 34016.134690, 35441.220090, 36884.850500, 38365.741930, 39903.026060, 41517.497110, 43230.364940, 45062.839420, 47040.785270, 49208.686640, 51615.682520, 54310.911900, 57345.705420, 60780.160230, 64676.565140, 69097.208930, 74081.734500, 79579.201070, 85516.021950, 91818.610450, 98416.457120, 105251.361400, 112268.200100, 119411.850000, 126629.198400, 133875.176600, 141106.726300, 148280.789300, 155377.514800, 162469.880800, 169654.072400, 177026.275100, 184730.628300, 193103.088700, 202527.567100, 213387.974400, 225960.352200, 240089.264300, 255511.405500, 271963.470200, 289091.769100, 306181.077400, 322425.786600, 337020.288000, 351006.047000, 372812.825500, 412717.459400, 480996.784700, 584604.351000, 717200.562800, 869122.538400, 1030707.396000};
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
    //loadmap[5] = 8372903.023;
    loadmap[5] = 14372903.023;
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
