#pragma once

#include "FEMSolveStrategy.h"
#include "../FEMDataType.h"
#include "../FEMMaterial.h"
#include "../FEMBoundary.h"
#include "../MatrixSolver/MatrixSolver.h"

#include <map>
#include <iostream>

using namespace std;

/*��������ͣ���ȡ��������������⡣����2d�������3d��ų��������ͨ��SolveStrategy���ı���ⷽʽ��NR/NDDR��
  ����Ӧ�ÿ�����չ��������
*/
class FEMSolver
{
public:
	FEMSolver();
	~FEMSolver();
	virtual void solve() = 0;

	void setSolveStrategy(FEMSolveStrategy* _strategy);
	void setMatrixSolver(MatrixSolver* const _matsolver);
	void setNodes(const int _numofnodes, CNode* const _nodes);
	void setVtxElements(const int _numofvtx, CVtxElement* const _vtxele);
	void setEdgElements(const int _numofedg, CEdgElement* const _edgele);
	void setTriElements(const int _numoftri, CTriElement* const _triele);
	void setMaterial(const std::map<int, FEMMaterial*> _materialmap);
	void setLoad(const std::map<int, double> _loadmap);
	void setBoundary(const std::map<int, FEMBoundary*> _boundarymap);
	void setMaxIterSteps(const int _maxitersteps);

protected:
	FEMSolveStrategy* strategy;
	MatrixSolver* matsolver;

	int m_num_nodes;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;

	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;

	int maxitersteps;
};

