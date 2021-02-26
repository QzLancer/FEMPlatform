#include "FEMMaterial.h"

#define PI 3.14159265358979323846

FEMMaterial::FEMMaterial()
{
    BHpoints = 0;
	mu = 4 * PI * 1e-7;
    Hdata = nullptr;
    Bdata = nullptr;
	linearflag = true;
}

FEMMaterial::~FEMMaterial()
{
	if (Hdata != NULL) {
		delete[] Hdata;
	}
	if (Bdata != NULL) {
		delete[] Bdata;
	}
}

double FEMMaterial::getMu(double B) const
{
    return mu;
}

double FEMMaterial::getdvdB(double B)
{
    double slope, H, b;
    if (BHpoints == 0 || B < 1e-9) return 0;
	getkHb(B, &slope, &H, &b);
    return -b / (B * B);
}

void FEMMaterial::getkHb(double B, double* k, double* H, double* b)
{
	if (B >= Bdata[BHpoints - 1]) {
		int  i = BHpoints - 2;
		(*k) = (Hdata[i] - Hdata[i - 1]) / (Bdata[i] - Bdata[i - 1]);
		(*b) = Hdata[i - 1] - (*k) * Bdata[i - 1];
	}
	else if (B < Bdata[0]) {
		(*k) = (Hdata[1] - Hdata[0]) / (Bdata[1] - Bdata[0]);
		(*b) = Hdata[0] - (*k) * Bdata[0];
	}
	else {
		for (int i = 0; i < BHpoints - 1; i++) {
			if (B >= Bdata[i] && B <= Bdata[i + 1]) {
				(*k) = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
				(*b) = Hdata[i] - (*k) * Bdata[i];
				break;
			}
		}
	}
	(*H) = (*k) * B + (*b);
}

void FEMMaterial::setName(const std::string name)
{
	this->name = name;
}

void FEMMaterial::setBHpoints(const int bhpoints)
{
	BHpoints = bhpoints;
}

void FEMMaterial::setBdata(double* const bdata)
{
	Bdata = new double[BHpoints];
	for (int i = 0; i < BHpoints; ++i) {
		Bdata[i] = bdata[i];
	}
}

void FEMMaterial::setHdata(double* const hdata)
{
	Hdata = hdata;
}

void FEMMaterial::setLinearFlag(const bool islinear)
{
	linearflag = islinear;
}

void FEMMaterial::setmu(const double mu)
{
	this->mu = mu;
}

bool FEMMaterial::getLinearFlag()
{
	return linearflag;
}

double* FEMMaterial::getBdata()
{
	return Bdata;
}

double* FEMMaterial::getHdata()
{
	return Hdata;
}
