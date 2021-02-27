#include "FEMModel.h"

FEMModel::DIMENSION FEMModel::getDimension() const
{
    return dimension;
}

std::string FEMModel::getMeshFile() const
{
    return meshfile;
}

std::vector<FEMMaterial*> FEMModel::getMaterialList() const
{
    return materiallist;
}

std::map<int, FEMMaterial*> FEMModel::getMaterialMap() const
{
    return materialmap;
}

std::map<int, double> FEMModel::getLoadMap() const
{
    return loadmap;
}

std::map<int, FEMBoundary*> FEMModel::getBoundaryMap() const
{
    return boundarymap;
}

FEMModel::FEMModel()
{

}

void FEMModel::init()
{
    setdimension();
    setMeshFile();
    createElement2Material();
    bulidGeometry2Load();
    buildGeometry2Constrain();
}
