#pragma once

#include "FEMModel.h"

//�������ŵĶ�ά�Ӵ���ģ�ͣ�����
class FEMContactorLinearModel :
	public FEMModel
{
protected:
	// ͨ�� FEMModel �̳�
	virtual void setModelName() override;
	virtual void setdimension() override;
	virtual void setMeshFile() override;
	virtual void createElement2Material() override;
	virtual void bulidGeometry2Load() override;
	virtual void buildGeometry2Constrain() override;
};

