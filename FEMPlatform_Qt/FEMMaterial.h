#pragma once

#define PI 3.14159265358979323846

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "FEMCoil.h"

#include <string>

/// @brief ���Բ��Ϲ������ṩFEM����Ҫ�Ĳ�������Ľӿڣ�����dv/dB�ȣ���֧��GPU����
class FEMMaterial
{
public:
	FEMMaterial();
	
	// ����������Ҫ�ͷ�CUDA�з�����ڴ�
	~FEMMaterial();

	/// @brief ͨ���ֶ����Բ�ֵ����ŵ���
	/// @param B �Ÿ�Ӧǿ��
	/// @return �ŵ���
	double getMu(double B = 0);

	/// @brief ͨ���ֶ����Բ�ֵ���������
	/// @param B �Ÿ�Ӧǿ��
	/// @return ������
	double getV(double B = 0);

	/// @brief ��������ʶԴŸ�Ӧǿ��ģ��ƫ��
	/// @param B �Ÿ�Ӧǿ��
	/// @return dvdB
	double getdvdB(double B);

	/// @brief ��������ʶԴŸ�Ӧǿ�ȵ�ƽ����ƫ��
	/// @param B �Ÿ�Ӧǿ��
	/// @return dvdB2 
	double getdvdB2(double B);

	/// @brief ͨ�����Բ�ֵ����BH����б��k���ų�ǿ��H�͸�б����B=0ʱ��Ӧ��H���洢Ϊb��
	/// @param B �Ÿ�Ӧǿ�ȣ��������
	/// @param k б�ʣ�H/B�����������
	/// @param H �ų�ǿ�ȣ��������
	/// @param b б��Ϊk���Ÿ�Ӧǿ��BΪ0ʱ��Ӧ�Ĵų�ǿ��
	void getkHb(double B, double* k, double* H, double* b);

	/// @brief 
	/// @return ������ 
	double getH_c() const;

	/// @brief 
	/// @return �������x��������ĽǶ�
	double getTheta_m() const;

	/// @brief 
	/// @return ��Ȧ����
	FEMCoil getFEMCoil() const;

	/// @brief ���ò�������
	/// @param name
	void setName(const std::string name);

	/// @brief ���÷����Բ��ϵ�BH����Ŀ
	/// @param bhpoints 
	void setBHpoints(const int bhpoints);

	/// @brief ���÷����Բ��ϵĴŸ�Ӧǿ�Ȳ�ֵ��
	/// @param bdata 
	void setBdata(double* const bdata);

	/// @brief ���÷����Բ��ϵĴų�ǿ�Ȳ�ֵ��
	/// @param hdata 
	void setHdata(double* const hdata);

	/// @brief ���ò����Ƿ�Ϊ����
	/// @param islinear 
	void setLinearFlag(const bool islinear);

	/// @brief �������Բ��ϵĴŵ���
	/// @param mu 
	void setMu(const double mu);

	/// @brief �������Ų��ϵĽ�����
	/// @param h_c 
	void setH_c(const double h_c);

	/// @brief �������Ų��ϵ�ʣ��
	/// @param b_r 
	void setB_r(const double b_r);

	/// @brief �������Ų��ϵĽǶ�
	/// @param theta_m 
	void setTheta_m(const double theta_m);

	/// @brief ������Ȧ����
	/// @param coil 
	void setFEMCoil(const FEMCoil coil);

	/// @brief 
	/// @return �жϲ����Ƿ�Ϊ���Բ���
	bool getLinearFlag();

	/// @brief 
	/// @return �����Բ��ϵĴŸ�Ӧǿ�Ȳ�ֵ��
	double* getBdata();

	/// @brief 
	/// @return �����Բ��ϵĴų�ǿ�Ȳ�ֵ��
	double* getHdata();

	/// @brief ͨ��GPU����ŵ��ʣ���ʵ�ַŵ�Դ�ļ��У������nvlink error����֪����ô������ʷ���ͷ�ļ���
	/// @param B �Ÿ�Ӧǿ��
	/// @return �ŵ���
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

	__device__ double getH_cinDevice()
	{
		return h_c;
	};

	__device__ double getTheta_minDevice()
	{
		return theta_m;
	};



	/// @brief ��CPU�д洢�Ĳ��ϲ���������GPU��
	/// @param material ��Ҫ�����Ĳ���
	void GPUCopy(FEMMaterial& material);

private:
	/// @brief ��������
	std::string name;
	/// @brief �ŵ���
	double mu;
	/// @brief ���ŵĽ��糡ǿ
	double h_c;	
	/// @brief ���ŵ�ʣ��
	double b_r;
	/// @brief �ų���x��������ļн�
	double theta_m;
	/// @brief �����Բ��ϵĴŸ�Ӧǿ������
	double* Bdata;
	/// @brief �����Բ��ϵĴų�ǿ������
	double* Hdata;
	/// @brief BH�����Ŀ
	int BHpoints;
	/// @brief �жϲ����Ƿ�Ϊ���Բ���
	bool linearflag;
	/// @brief �ж�Bdata��Hdata�Ƿ����ΪUnified Memory����ѡ��������ʽ
	bool gpuflag;
	/// @brief ��Ȧ����
	FEMCoil coil;
};