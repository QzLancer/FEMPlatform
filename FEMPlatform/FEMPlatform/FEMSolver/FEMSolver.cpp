#include "FEMSolver.h"

FEMSolver::FEMSolver():
    strategy(nullptr),
    matsolver(nullptr),
    m_num_nodes(0),
    m_num_vtxele(0),
    m_num_edgele(0),
    m_num_triele(0),
    mp_node(nullptr),
    mp_vtxele(nullptr),
    mp_edgele(nullptr),
    mp_triele(nullptr)
{

}

FEMSolver::~FEMSolver()
{
    if (mp_triele != nullptr)
        delete[] mp_triele;
    if (mp_edgele != nullptr)
        delete[] mp_edgele;
    if (mp_vtxele != nullptr)
        delete[] mp_vtxele;
    if (mp_node != nullptr)
        delete[] mp_node;
    if(strategy != nullptr)
        delete strategy;
}

void FEMSolver::setSolveStrategy(FEMSolveStrategy* _strategy)
{
    strategy = _strategy;
}

void FEMSolver::setMatrixSolver(MatrixSolver* const _matsolver)
{
    matsolver = _matsolver;
}

void FEMSolver::setNodes(const int _numofnodes, CNode* const _nodes)
{
    m_num_nodes = _numofnodes;
    mp_node = _nodes;
}

void FEMSolver::setVtxElements(const int _numofvtx, CVtxElement* const _vtxele)
{
    m_num_vtxele = _numofvtx;
    mp_vtxele = _vtxele;
}

void FEMSolver::setEdgElements(const int _numofedg, CEdgElement* const _edgele)
{
    m_num_edgele = _numofedg;
    mp_edgele = _edgele;
}

void FEMSolver::setTriElements(const int _numoftri, CTriElement* const _triele)
{
    m_num_triele = _numoftri;
    mp_triele = _triele;
}

void FEMSolver::setMaterial(const std::map<int, FEMMaterial*> _materialmap)
{
    materialmap = _materialmap;
}

void FEMSolver::setLoad(const std::map<int, double> _loadmap)
{
    loadmap = _loadmap;
}

void FEMSolver::setBoundary(const std::map<int, FEMBoundary*> _boundarymap)
{
    boundarymap = _boundarymap;
}

void FEMSolver::setMaxIterSteps(const int _maxitersteps)
{
    maxitersteps = _maxitersteps;
}
