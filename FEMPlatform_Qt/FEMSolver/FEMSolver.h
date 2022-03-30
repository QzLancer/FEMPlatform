 #pragma once

#include "FEMSolveStrategy.h"
#include "../FEMDataType.h"
#include "../FEMMaterial.h"
#include "../FEMBoundary.h"
#include "../MatrixSolver/MatrixSolver.h"
#include "../FEMModel/FEMModel.h"
#include "MeshManager/FEMMeshManager.h"

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <math.h>
#include "..\FEMMovingPart.h"

#define PI 3.14159265358979323846

using namespace std;

/*��������ͣ���ȡ��������������⡣����2d�������3d��ų��������ͨ��SolveStrategy���ı���ⷽʽ��NR/NDDR��
  ����Ӧ�ÿ�����չ��������
*/
class FEMSolver
{
public:
	FEMSolver();
	virtual ~FEMSolver();
	virtual void solveStatic() = 0;
	virtual void solveDynamic() = 0;
	virtual void solveMagneticForce() = 0;
	virtual void solveMagneticForce1() = 0;

	void setSolveStrategy(FEMSolveStrategy* _strategy);
	void setMatrixSolver(MatrixSolver* const _matsolver);
	void setModelName(const std::string _modelname);
	virtual void setNodes(const int _numofnodes, CNode* const _nodes);
	virtual void setVtxElements(const int _numofvtx, CVtxElement* const _vtxele);
	virtual void setEdgElements(const int _numofedg, CEdgElement* const _edgele);
	virtual void setTriElements(const int _numoftri, CTriElement* const _triele);
	void setMaterial(const std::map<int, FEMMaterial*> _materialmap);
	void setLoad(const std::map<int, double> _loadmap);
	void setBoundary(const std::map<int, FEMBoundary*> _boundarymap);
	void setMaxIterSteps(const int _maxitersteps);
	void setMaxError(const double _error);
	void writeVtkFile(std::string _name);
	void writeTxtFile(std::string _name);
	void writeVtkFileNoAir(std::string _name, vector<int> air_domain);
	void writeGeometryVtkFile(std::string _name);
	virtual void setDeformedDomain(const std::vector<int> _deformedlist);
	virtual void setMovingPart(const std::map<int, FEMMovingPart*> _movingmap);
	void setMeshManager(FEMMeshManager* _meshmanager);
	void updateLoadmap(int domain, double current);

	//std::vector<double> getA() const;
	FEMModel::DIMENSION dimension;

protected:
	virtual void processBoundaryCondition() = 0;
	virtual void processMaterial() = 0;
	virtual void processLoad() = 0;

	FEMSolveStrategy* strategy;
	MatrixSolver* matsolver;

	int m_num_nodes;
	int m_num_vtxele;
	int m_num_edgele;
	int m_num_triele;

	std::string modelname;
	CNode* mp_node;
	CVtxElement* mp_vtxele;
	CEdgElement* mp_edgele;
	CTriElement* mp_triele;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
	std::vector<int> deformedlist;	//�α�����
	std::map<int, FEMMovingPart*>movingmap;	//�˶�����

	int maxitersteps;
	double maxerror;

	int num_freenodes;	//���ɽڵ���Ŀ
	std::vector<int> node_reorder;	//ǰnum_dof��Ԫ�ض�Ӧ�Ǳ߽�ڵ㣬֮���Ԫ�ض�Ӧ��һ��߽�����
	std::vector<int> node_pos;	//ԭ�ڵ��Ŷ�Ӧ��reorder��Ľڵ���

	//std::vector<double> A{ 0 };
	//std::vector<double> At{ 0 };
	//std::vector<double> Bx{ 0 }, By{ 0 }, Bz{ 0 }, B{ 0 };
	FEMMeshManager* meshmanager;
	double* current;
	double* dis;
	double* velocity;
	double* acc;
	double* magneticforce;
	double* springforce;
	double* flux;
	double mass;
	double h, U, R;
	int staticsteps;
	int dynamicsteps;
};

