#include "FEMMaterial.h"

FEMMaterial::FEMMaterial():
	BHpoints(0),
	mu(PI * 4e-7),
	Hdata(nullptr),
	Bdata(nullptr),
	linearflag(true)
{
}

FEMMaterial::~FEMMaterial()
{
	if (gpuflag == true) {
		if (Hdata != nullptr) {
			cudaFree(Hdata);
			Hdata = nullptr;
		}
		if (Bdata != nullptr) {
			cudaFree(Bdata);
			Bdata = nullptr;
		}
	}
	else {
		if (Hdata != nullptr) {
			delete[] Hdata;
			Hdata = nullptr;
		}
		if (Bdata != nullptr) {
			delete[] Bdata;
			Bdata = nullptr;
		}
	}
}

double FEMMaterial::getMu(double B)
{
	//double slope, H, b;

	//if (linearflag == false) {
	//	if (B < 1e-3) {
	//		mu = Bdata[1] / Hdata[1];
	//	}
	//	else {
	//		getkHb(B, &slope, &H, &b);
	//		if (B / H < PI * 4e-7) {
	//			mu = PI * 4e-7;
	//		}
	//		else {
	//			mu = B / H;
	//		}
	//	}
	//}
	//return mu;

	double slope, H, b;
	if (linearflag == true) {
		return mu;
	}
	else {
		if (B < 1e-3)  return Bdata[1] / Hdata[1];
		getkHb(B, &slope, &H, &b);
	}
	if (B / H < PI * 4e-7)  return PI * 4e-7;
	return B / H;
}

double FEMMaterial::getdvdB(double B)
{
    double slope, H, b;
    if (BHpoints == 0 || B < 1e-9) return 0;
	getkHb(B, &slope, &H, &b);
    return -b / (B * B);
}

//线性插值计算斜率k，磁场强度H
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
void FEMMaterial::GPUCopy(FEMMaterial& material)
{
	name = material.name;
	mu = material.mu;
	BHpoints = material.BHpoints;
	linearflag = material.linearflag;
	gpuflag = true;
	if (linearflag == false) {
		cudaMallocManaged((void**)&Bdata, BHpoints * sizeof(double));
		memcpy(Bdata, material.Bdata, BHpoints * sizeof(double));
		cudaMallocManaged((void**)&Hdata, BHpoints * sizeof(double));
		memcpy(Hdata, material.Hdata, BHpoints * sizeof(double));
	}
}
//
//__device__ inline double FEMMaterial::getMuinDevice(double B)
//{
//	return 0;
//}
//
//__device__ inline double FEMMaterial::getdvdBinDevice(double B)
//{
//	return 0;
//}
//
//__device__ inline void FEMMaterial::getkHbinDevice(double B, double* k, double* H, double* b)
//{
//	
//}
//
//__device__ inline bool FEMMaterial::getLinearFlaginDevice()
//{
//	return true;
//}
