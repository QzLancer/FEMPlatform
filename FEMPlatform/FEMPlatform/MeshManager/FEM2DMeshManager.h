#pragma once
#include "FEMMeshManager.h"

#include <vector>

//Ŀǰ����ģ�͵Ĵ���������������ͨ��ģ�����ݴ���͵���ģ�͵ı߽硢����������ӷ���
class FEM2DGenerator :
    public FEMMeshManager
{
public:
    ~FEM2DGenerator();
    virtual void readMeshFile(string meshfile) override;

private:
    void read2DMphtxt(string meshfile); //��ȡCOMSOL5.5�汾�����ļ�
};

