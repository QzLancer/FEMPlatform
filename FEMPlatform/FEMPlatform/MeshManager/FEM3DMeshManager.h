#pragma once
#include "FEMMeshManager.h"
class FEM3DMeshManager :
    public FEMMeshManager
{
	virtual void readMeshFile(string meshfile) override;
};

