#pragma once
#include "FEMGenerator.h"

#include <vector>

//目前单个模型的处理放在这里，后续将通用模型数据处理和单独模型的边界、负载条件添加分离
class FEM2DGenerator :
    public FEMGenerator
{
public:
    ~FEM2DGenerator();
    virtual void readMeshFile(string meshfile) override;
    virtual void createElement2Material() override;
    virtual void bulidGeometry2Load() override;
    virtual void buildGeometry2Constrain() override;

private:
    void read2DMphtxt(string meshfile);
};

