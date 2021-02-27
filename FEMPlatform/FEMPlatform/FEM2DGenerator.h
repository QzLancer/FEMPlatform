#pragma once
#include "FEMGenerator.h"

#include <vector>

//Ŀǰ����ģ�͵Ĵ���������������ͨ��ģ�����ݴ���͵���ģ�͵ı߽硢����������ӷ���
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

