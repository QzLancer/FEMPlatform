#pragma once

#include "FEMModel.h"

class FEMContactorNonLinearModel :
	public FEMModel
{
	// Í¨¹ý FEMModel ¼Ì³Ð
	virtual void setModelName() override;
	virtual void setdimension() override;
	virtual void setMeshFile() override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
};

