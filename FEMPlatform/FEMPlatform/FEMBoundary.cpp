#include "FEMBoundary.h"

std::string FEMAxialSymmetry::getName() const
{
    return "AxialSymmetry";
}

std::string FEMMagneticInsulation::getName() const
{
    return "MagneticInsulation";
}
