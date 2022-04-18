#include "RelayModelwithoutBand.h"

void RelayModelwithoutBand::setModelName()
{
	modelname = "RelayModelwithoutBand";
}

void RelayModelwithoutBand::setDimension()
{
	dimension = DIMENSION::D2AXISM;
}

void RelayModelwithoutBand::setFile()
{
	//geofile = "D:/femplatform/model/qzlancer_model2.geo";
	meshfile = "../model/modelwithband_0.mphtxt";
}

void RelayModelwithoutBand::createElement2Material()
{
	addLinearMaterial("Air", 4 * PI * 1e-7);
	addLinearMaterial("Iron", 300 * 4 * PI * 1e-7);
	double* bdata = new double[] { 0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
	double* hdata = new double[] { 0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
	addNonlinearMaterial("Soft Iron1", 16, bdata, hdata);

	FEMCoil coil;
	coil.xyz = 1;
	coil.direction = 1;
	coil.Nc = 2400;
	coil.Rc = 40;
	coil.Vc = 24;
	coil.type = CICULARCOIL;
	coil.height = 0.0147;
	coil.center[0] = 0;
	coil.center[1] = 0;
	coil.center[2] = 0;
	coil.r0 = 0;
	coil.width = 0.011687;
	coil.Jor = coil.Vc / coil.Rc * coil.Nc / coil.width / coil.height;
	coil.tau = coil.Nc / coil.width / coil.height;
	addCoil("Coil", coil);

	//materialmap[1] = materiallist[2];   //ÏÎÌú
	//materialmap[2] = materiallist[2];   //éîÌú
	//materialmap[3] = materiallist[3];   //ÏßÈ¦
	//materialmap[4] = materiallist[0];   //¿ÕÆø
	//materialmap[5] = materiallist[0];   //Ñ¹Ëõ¿ÕÆø

	materialmap[1] = materiallist[0];   //¿ÕÆø
	materialmap[2] = materiallist[0];   //¿ÕÆø
	materialmap[3] = materiallist[0];   //Ñ¹Ëõ¿ÕÆø
	materialmap[4] = materiallist[2];   //ÏÎÌú
	materialmap[5] = materiallist[2];   //éîÌú
	materialmap[6] = materiallist[0];   //¿ÕÆø
	materialmap[7] = materiallist[3];   //ÏßÈ¦
}

void RelayModelwithoutBand::bulidGeometry2Load()
{
	//loadmap[3] = materialmap[3]->getFEMCoil().Jor;
	loadmap[7] = materialmap[7]->getFEMCoil().Jor;
}

void RelayModelwithoutBand::buildGeometry2Constrain()
{
	//std::vector<int> id = { 33, 34, 35, 36};
	//for (auto a : id) {
	//	boundarymap[a] = new FEMAxialSymmetry;
	//}
	//std::vector<int> id1 = { 31, 32 };
	//for (auto b : id1) {
	//	boundarymap[b] = new FEMMagneticInsulation;
	//}

	std::vector<int> id = { 1, 2, 3, 5, 6, 8 };
	for (auto a : id) {
		boundarymap[a] = new FEMAxialSymmetry;
	}
	std::vector<int> id1 = { 40, 43 };
	for (auto b : id1) {
		boundarymap[b] = new FEMMagneticInsulation;
	}
}

void RelayModelwithoutBand::setUnitRatio()
{
	unitratio = 1;
}

void RelayModelwithoutBand::buildGeometry2Deformed()
{
	deformedlist.push_back(3);
}

void RelayModelwithoutBand::buildGeometry2MovingPart()
{
	FEMMovingPart* moving1 = new FEMMovingPart;
	moving1->direction[1] = 1;
	moving1->limit[1].max = 0.00249;
	moving1->setMass(0.024);
	double* pos = new double[4]{ 0, 0.0016999, 0.0017, 0.0027 };
	double* force = new double[4]{ -6.0, -6.63, -13.63, -27.0 };
	moving1->setSpringForce(4, pos, force);

	movingmap[4] = moving1;
}
