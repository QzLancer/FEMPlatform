#pragma once

#include "FEMModel.h"

class FEMContactorNonLinearModel :
	public FEMModel
{
	// ͨ�� FEMModel �̳�
	virtual void setModelName() override;
	virtual void setDimension() override;
	virtual void setFile() override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
	virtual void setUnitRatio() override;
};

