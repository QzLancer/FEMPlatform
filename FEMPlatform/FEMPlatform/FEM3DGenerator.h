#pragma once
#include "FEMGenerator.h"
class FEM3DGenerator :
    public FEMGenerator
{
	virtual void readMeshFile(string meshfile) override;
	virtual void createElement2Material() override;	/*设置单元材料类型*/
	virtual void bulidGeometry2Load() override;		/*设置负载*/
	virtual void buildGeometry2Constrain() override;	/*设置边界条件*/
};

