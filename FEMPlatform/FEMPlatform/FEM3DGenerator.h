#pragma once
#include "FEMGenerator.h"
class FEM3DGenerator :
    public FEMGenerator
{
	virtual void readMeshFile(string meshfile) override;
	virtual void createElement2Material() override;	/*���õ�Ԫ��������*/
	virtual void bulidGeometry2Load() override;		/*���ø���*/
	virtual void buildGeometry2Constrain() override;	/*���ñ߽�����*/
};

