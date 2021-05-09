#pragma once
#include "FEMMeshManager.h"

#include <vector>

//Ŀǰ����ģ�͵Ĵ���������������ͨ��ģ�����ݴ���͵���ģ�͵ı߽硢����������ӷ���
class FEM2DMeshManager :
    public FEMMeshManager
{
public:
    ~FEM2DMeshManager();
    virtual void readMeshFile(string meshfile) override;

private:
    void read2DMphtxt(string meshfile); //��ȡCOMSOL5.5�汾�����ļ�
    void read2DMsh(string meshfile = "");    //��ȡGmsh�����ļ�
};

