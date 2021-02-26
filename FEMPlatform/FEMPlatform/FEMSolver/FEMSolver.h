#pragma once

#include "FEMSolveStrategy.h"

class FEMSolver
{
public:
	void setSolveStrategy(FEMSolveStrategy* strategy);

protected:
	FEMSolveStrategy* strategy;
};

