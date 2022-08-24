#pragma once

#define PI 3.14159265358979323846

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "FEMCoil.h"

#include <string>

/// @brief 磁性材料管理，并提供FEM所需要的参数计算的接口（例如dv/dB等），支持GPU计算
class FEMMaterial
{
public:
	FEMMaterial();
	
	// 析构函数需要释放CUDA中分配的内存
	~FEMMaterial();

	/// @brief 通过分段线性插值计算磁导率
	/// @param B 磁感应强度
	/// @return 磁导率
	double getMu(double B = 0);

	/// @brief 通过分段线性插值计算磁阻率
	/// @param B 磁感应强度
	/// @return 磁阻率
	double getV(double B = 0);

	/// @brief 计算磁阻率对磁感应强度模的偏导
	/// @param B 磁感应强度
	/// @return dvdB
	double getdvdB(double B);

	/// @brief 计算磁阻率对磁感应强度的平方的偏导
	/// @param B 磁感应强度
	/// @return dvdB2 
	double getdvdB2(double B);

	/// @brief 通过线性插值计算BH曲线斜率k、磁场强度H和该斜率下B=0时对应的H（存储为b）
	/// @param B 磁感应强度，输入参数
	/// @param k 斜率（H/B），输出参数
	/// @param H 磁场强度，输出参数
	/// @param b 斜率为k，磁感应强度B为0时对应的磁场强度
	void getkHb(double B, double* k, double* H, double* b);

	/// @brief 
	/// @return 矫顽力 
	double getH_c() const;

	/// @brief 
	/// @return 永磁相对x轴正半轴的角度
	double getTheta_m() const;

	/// @brief 
	/// @return 线圈参数
	FEMCoil getFEMCoil() const;

	/// @brief 设置材料名称
	/// @param name
	void setName(const std::string name);

	/// @brief 设置非线性材料的BH点数目
	/// @param bhpoints 
	void setBHpoints(const int bhpoints);

	/// @brief 设置非线性材料的磁感应强度插值点
	/// @param bdata 
	void setBdata(double* const bdata);

	/// @brief 设置非线性材料的磁场强度插值点
	/// @param hdata 
	void setHdata(double* const hdata);

	/// @brief 设置材料是否为线性
	/// @param islinear 
	void setLinearFlag(const bool islinear);

	/// @brief 设置线性材料的磁导率
	/// @param mu 
	void setMu(const double mu);

	/// @brief 设置永磁材料的矫顽力
	/// @param h_c 
	void setH_c(const double h_c);

	/// @brief 设置永磁材料的剩磁
	/// @param b_r 
	void setB_r(const double b_r);

	/// @brief 设置永磁材料的角度
	/// @param theta_m 
	void setTheta_m(const double theta_m);

	/// @brief 设置线圈属性
	/// @param coil 
	void setFEMCoil(const FEMCoil coil);

	/// @brief 
	/// @return 判断材料是否为线性材料
	bool getLinearFlag();

	/// @brief 
	/// @return 非线性材料的磁感应强度插值点
	double* getBdata();

	/// @brief 
	/// @return 非线性材料的磁场强度插值点
	double* getHdata();

	/// @brief 通过GPU计算磁导率，把实现放到源文件中，会出现nvlink error，不知道怎么解决，故放在头文件中
	/// @param B 磁感应强度
	/// @return 磁导率
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



	/// @brief 将CPU中存储的材料参数拷贝到GPU中
	/// @param material 需要拷贝的材料
	void GPUCopy(FEMMaterial& material);

private:
	/// @brief 材料名称
	std::string name;
	/// @brief 磁导率
	double mu;
	/// @brief 永磁的矫顽场强
	double h_c;	
	/// @brief 永磁的剩磁
	double b_r;
	/// @brief 磁场与x轴正方向的夹角
	double theta_m;
	/// @brief 非线性材料的磁感应强度数组
	double* Bdata;
	/// @brief 非线性材料的磁场强度数组
	double* Hdata;
	/// @brief BH点的数目
	int BHpoints;
	/// @brief 判断材料是否为线性材料
	bool linearflag;
	/// @brief 判断Bdata和Hdata是否分配为Unified Memory，以选择析构方式
	bool gpuflag;
	/// @brief 线圈属性
	FEMCoil coil;
};