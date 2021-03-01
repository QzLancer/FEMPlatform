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

std::string FEMModel::getModelName() const
{
    return modelname;
}

FEMModel::FEMModel()
{

}

void FEMModel::init()
{
    setModelName();
    setdimension();
    setMeshFile();
    createElement2Material();
    bulidGeometry2Load();
    buildGeometry2Constrain();
}

void FEMModel::addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata)
{
    FEMMaterial* material = new FEMMaterial;
    material->setName(_name);
    material->setLinearFlag(false);
    material->setBHpoints(_bhpoints);
    material->setBdata(_bdata);
    material->setHdata(_hdata);
    materiallist.push_back(material);
}

void FEMModel::addLinearMaterial(std::string _name, double _mu)
{
    FEMMaterial* material = new FEMMaterial;
    material->setBHpoints(1);
    material->setLinearFlag(true);
    material->setmu(_mu);
    materiallist.push_back(material);
}
