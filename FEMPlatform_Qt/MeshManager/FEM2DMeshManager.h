#pragma once
#include "FEMMeshManager.h"

#include <vector>

//目前单个模型的处理放在这里，后续将通用模型数据处理和单独模型的边界、负载条件添加分离
class FEM2DMeshManager :
    public FEMMeshManager
{
public:
    ~FEM2DMeshManager();
    virtual void readMeshFile(string meshfile) override;

private:
    void read2DMphtxt(string meshfile); //读取COMSOL5.5版本分网文件
    void read2DMsh(string meshfile = "");    //读取Gmsh分网文件
};

