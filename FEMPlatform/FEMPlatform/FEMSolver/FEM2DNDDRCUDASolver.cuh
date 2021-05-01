#pragma once
#include "FEM2DSolver.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class FEM2DNDDRCUDASolver :
    public FEM2DSolver
{
public:
    virtual void solveStatic() override;
	virtual ~FEM2DNDDRCUDASolver() override;

	//virtual void setNodes(const int _numofnodes, CNode* const _nodes) override;
	//virtual void setVtxElements(const int _numofvtx, CVtxElement* const _vtxele) override;
	//virtual void setEdgElements(const int _numofedg, CEdgElement* const _edgele) override;
	//virtual void setTriElements(const int _numoftri, CTriElement* const _triele) override;
	virtual void processMaterial() override;    //将材料添加到单元上

	void GPUFree();

protected:

private:
	int CudaThrdNum = 300;
	int CudaBlckNum = 256;

	CNode* d_mp_node;
	CTriElement* d_mp_triele;
	FEMMaterial* materialarray;
	FEMMaterial* d_materialarray;

	void GPUInitialMallocCopy();
	void processNDDRNode();
};

__global__ void nodeAnalysis(int d_m_num_nodes, CNode* d_mp_node, CTriElement* d_mp_triele);

__global__ void calculateGlobalError(int d_m_num_nodes, CNode* d_mp_node, CTriElement* d_mp_triele, double* a, double* b);

__global__ void copyAttoAtold(int d_m_num_nodes, CNode* d_mp_node);

__global__ void makeTrangleinDevice(int numofTrangle, CNode* d_mp_node, CTriElement* d_mp_triele);

__global__ void assignMattoTriEle(int numofTrangle, CTriElement* d_mp_triele, FEMMaterial* d_materialarray);
