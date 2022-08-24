#pragma once

#include <string>

/// @brief 边界条件类，目前只支持轴对称边界和自由边界条件（有限元分析中采用将节点自由度直接赋值为0的计算方法）
class FEMBoundary
{
public:
	virtual std::string getName() const = 0;
	int getBoundaryType() const;

protected:
	int boundarytype;	//边界条件，用于确认自由度

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