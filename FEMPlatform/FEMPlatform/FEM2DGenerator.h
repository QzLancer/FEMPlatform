#pragma once
#include "FEMGenerator.h"

#include <vector>
class FEM2DGenerator :
    public FEMGenerator
{
public:
    virtual void readMeshFile(string meshfile) override;
    virtual void createElement2Material() override;
    virtual void bulidGeometry2Load() override;
    virtual void buildGeometry2Constrain() override;

private:
    int m_num_nodes;
    int m_num_vtxele;
    int m_num_edgele;
    int m_num_triele;
    C2DNode* mp_2Dnode;
    CVtxElement* mp_vtxele;
    CEdgElement* mp_edgele;
    CTriElement* mp_triele;

    void read2DMphtxt(string meshfile);
};

