#pragma once

#include <string>

class FEMBoundary
{
public:
	virtual std::string getName() const = 0;

};

class FEMAxialSymmetry :
	public FEMBoundary 
{
public:
	virtual std::string getName() const override;
};

class FEMMagneticInsulation :
	public FEMBoundary
{
public:
	virtual std::string getName() const override;
};