#include "FEMBoundary.h"

FEMAxialSymmetry::FEMAxialSymmetry()
{
    boundarytype = 1;
}

std::string FEMAxialSymmetry::getName() const
{
    return "AxialSymmetry";
}

FEMMagneticInsulation::FEMMagneticInsulation()
{
    boundarytype = 1;
}

std::string FEMMagneticInsulation::getName() const
{
    return "MagneticInsulation";
}

int FEMBoundary::getBoundaryType() const
{
    return boundarytype;
}
