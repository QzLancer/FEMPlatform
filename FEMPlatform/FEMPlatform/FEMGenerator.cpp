#include "FEMGenerator.h"


void FEMGenerator::addNonlinearMaterial(string name, int bhpoints, double* bdata, double* hdata)
{
	FEMMaterial* material = new FEMMaterial;
	material->setName(name);
	material->setLinearFlag(false);
	material->setBHpoints(bhpoints);
	material->setBdata(bdata);
	material->setHdata(hdata);
	materiallist.push_back(material);
}

void FEMGenerator::addLinearMaterial(string name, double mu)
{
	FEMMaterial* material = new FEMMaterial;
	material->setBHpoints(1);
	material->setLinearFlag(true);
	material->setmu(mu);
}
