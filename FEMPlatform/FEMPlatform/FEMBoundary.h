#pragma once

#include <string>

class FEMBoundary
{
public:
	virtual std::string getName() const = 0;
	int getBoundaryType() const;

protected:
	int boundarytype;	//�߽�����������ȷ�����ɶ�

};

class FEMAxialSymmetry :
	public FEMBoundary 
{
public:
	FEMAxialSymmetry();
	virtual std::string getName() const override;
};

class FEMMagneticInsulation :
	public FEMBoundary
{
public:
	FEMMagneticInsulation();
	virtual std::string getName() const override;
};