#pragma once

#define PI 3.14159265358979323846

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <string>

class FEMMaterial
{
public:
	FEMMaterial();
	~FEMMaterial();
	double getMu(double B = 0);
	double getdvdB(double B) ;	//���Դ���
	void getkHb(double B, double* k, double* H, double* b);
	void setName(const std::string name);
	void setBHpoints(const int bhpoints);
	void setBdata(double* const bdata);
	void setHdata(double* const hdata);
	void setLinearFlag(const bool islinear);
	void setmu(const double mu);

	bool getLinearFlag();
	double* getBdata();
	double* getHdata();
	
	//��ʵ�ַŵ�Դ�ļ��У������nvlink error����֪����ô���
	__device__ double getMuinDevice(double B = 0)
	{
		double slope, H, b;
		//printf("__device__ double getMuinDevice, B = %f\n", B);
		if (linearflag == true) {
			return mu;
		}
		else {
			if (B < 1e-3)  return Bdata[1] / Hdata[1];
			getkHbinDevice(B, &slope, &H, &b);
		}
		if (B / H < PI * 4e-7)  return PI * 4e-7;
		return B / H;
	};

	__device__ void getkHbinDevice(double B, double* k, double* H, double* b)
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
	};

	__device__ double getdvdBinDevice(double B)
	{
		double slope, H, b;
		if (BHpoints == 0 || B < 1e-9) return 0;
		getkHbinDevice(B, &slope, &H, &b);
		return -b / (B * B);
	};


	__device__ bool getLinearFlaginDevice()
	{
		return linearflag;
	};

	//GPU
	void GPUCopy(FEMMaterial& material);
private:
	std::string name;
	double mu;
	double* Bdata;
	double* Hdata;
	int BHpoints;
	bool linearflag;
	bool gpuflag{ false };	//�ж�Bdata��Hdata�Ƿ����ΪUnified Memory����ѡ��������ʽ
};

