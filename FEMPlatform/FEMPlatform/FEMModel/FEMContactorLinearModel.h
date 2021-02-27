#pragma once

#include "FEMModel.h"

//�������ŵĶ�ά�Ӵ���ģ�ͣ�����
class FEMContactorLinearModel :
	public FEMModel
{
protected:
	// ͨ�� FEMModel �̳�
	virtual void setdimension() override;
	virtual void setMeshFile() override;
	virtual void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata) override;
	virtual void addLinearMaterial(std::string _name, double _mu) override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
};
