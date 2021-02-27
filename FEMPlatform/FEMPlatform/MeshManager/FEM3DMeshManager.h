#pragma once
#include "FEMMeshManager.h"
class FEM3DGenerator :
    public FEMMeshManager
{
	virtual void readMeshFile(string meshfile) override;
};

