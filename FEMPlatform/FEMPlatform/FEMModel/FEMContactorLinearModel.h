#pragma once

#include "FEMModel.h"

//不含永磁的二维接触器模型，线性
class FEMContactorLinearModel :
	public FEMModel
{
protected:
	// 通过 FEMModel 继承
	virtual void setModelName() override;
	virtual void setDimension() override;
	virtual void setFile() override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
	virtual void setUnitRatio() override;
};

