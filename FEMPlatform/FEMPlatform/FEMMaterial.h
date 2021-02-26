#pragma once
#include <string>

class FEMMaterial
{
public:
	FEMMaterial();
	~FEMMaterial();
	double getMu(double B) const;
	double getdvdB(double B) ;
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
	
private:
	std::string name;
	double mu;
	double* Bdata;
	double* Hdata;
	int BHpoints;
	bool linearflag;
};

