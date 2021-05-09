#include "RelayDynamicModel.h"

void RelayDynamicModel::setModelName()
{
	modelname = "RelayDynamic";
}

void RelayDynamicModel::setDimension()
{
	dimension = DIMENSION::D2AXISM;
}

void RelayDynamicModel::setFile()
{
	geofile = "D:/femplatform/model/modelwithband.geo";
}

void RelayDynamicModel::createElement2Material()
{
	addLinearMaterial("Air", 4 * PI * 1e-7);
	double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
	double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
	addNonlinearMaterial("Soft Iron1", 16, bdata, hdata);

	materialmap[1] = materiallist[1];   //ÏÎÌú
	materialmap[2] = materiallist[1];   //éîÌú
	materialmap[3] = materiallist[0];   //ÏßÈ¦
	materialmap[4] = materiallist[0];   //¿ÕÆø
	materialmap[5] = materiallist[0];   //Ñ¹Ëõ¿ÕÆø
	materialmap[6] = materiallist[0];   //¿ÕÆø
	materialmap[7] = materiallist[0];   //¿ÕÆø
}

void RelayDynamicModel::bulidGeometry2Load()
{
	loadmap[3] = 8e6;
}

void RelayDynamicModel::buildGeometry2Constrain()
{
	std::vector<int> id = { 33, 34, 35, 36, 41, 42 };
	for (auto a : id) {
		boundarymap[a] = new FEMAxialSymmetry;
	}
	std::vector<int> id1 = { 31, 32 };
	for (auto b : id1) {
		boundarymap[b] = new FEMMagneticInsulation;
	}
}

void RelayDynamicModel::setUnitRatio()
{
	unitratio = 1;
}

void RelayDynamicModel::buildGeometry2Deformed()
{
	deformedlist.push_back(5);
}

void RelayDynamicModel::buildGeometry2MovingPart()
{
	FEMMovingPart moving1;
	moving1.direction[1] = 1;
	moving1.limit[1].max = 0.002;
	movingmap[1] = moving1;
}
