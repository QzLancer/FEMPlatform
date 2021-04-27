#include "FEMModel.h"

#include <iostream>

FEMModel::DIMENSION FEMModel::getDimension() const
{
    return dimension;
}

std::string FEMModel::getMeshFile() const
{
    return meshfile;
}

std::string FEMModel::getGeoFile() const
{
    return geofile;
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

double FEMModel::getUnitRatio() const
{
    return unitratio;
}

FEMModel::FEMModel()
{

}

void FEMModel::init()
{
    setModelName();
    setDimension();
    setFile();
    createElement2Material();
    bulidGeometry2Load();
    buildGeometry2Constrain();
    setUnitRatio();
    buildGeometry2Deformed();
    buildGeometry2Constrain();
    buildGeometry2MovingPart();

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
    material->setName(_name);
    material->setBHpoints(0);
    material->setLinearFlag(true);
    material->setMu(_mu);
    materiallist.push_back(material);
}

void FEMModel::addPermMagentMaterial(std::string _name, double _mu, double _h_c, double _b_r, double _theta_m)
{
    FEMMaterial* material = new FEMMaterial;
    material->setName(_name);
    material->setLinearFlag(true);
    material->setMu(_mu);
    material->setBHpoints(0);
    material->setH_c(_h_c);
    material->setB_r(_b_r);
    material->setTheta_m(_theta_m);
    materiallist.push_back(material);
}

void FEMModel::addCoil(std::string _name, FEMCoil _coil)
{
    FEMMaterial* material = new FEMMaterial;
    material->setName(_name);
    material->setBHpoints(0);
    material->setLinearFlag(true);
    material->setMu(4 * PI * 1e-7);
    material->setFEMCoil(_coil);
    materiallist.push_back(material);
}
