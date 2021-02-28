#pragma once

#include "FEMModel.h"

//不含永磁的二维接触器模型，线性
class FEMContactorLinearModel :
	public FEMModel
{
protected:
	// 通过 FEMModel 继承
	virtual void setModelName() override;
	virtual void setdimension() override;
	virtual void setMeshFile() override;
	virtual void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata) override;
	virtual void addLinearMaterial(std::string _name, double _mu) override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
};

