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

/*求解器类型，读取求解参数，进行求解。派生2d求解器和3d电磁场求解器，通过SolveStrategy来改变求解方式（NR/NDDR）
  后续应该考虑拓展其他物理场
*/
class FEMSolver
{
public:
	FEMSolver();
	virtual ~FEMSolver();

	/// @brief 静态特性计算
	virtual void solveStatic() = 0;

	/// @brief 动态特性计算
	virtual void solveDynamic() = 0;

	/// @brief 计算电磁力
	virtual void solveMagneticForce() = 0;
	virtual void solveMagneticForce1() = 0;

	/// @brief 设置求解策略，目前没用上
	/// @param _strategy 
	void setSolveStrategy(FEMSolveStrategy* _strategy);

	/// @brief 设置矩阵求解器
	/// @param _matsolver 
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

	/// @brief 求解结果写入vtk文件
	/// @param _name 文件名称
	void writeVtkFile(std::string _name);

	/// @brief 求解结果写入txt文件
	/// @param _name 文件名称
	void writeTxtFile(std::string _name);

	/// @brief 除空气域之外的求解结果写入vtk文件
	/// @param _name 文件名称
	/// @param air_domain 空气域
	void writeVtkFileNoAir(std::string _name, vector<int> air_domain);

	/// @brief 将几何写入vtk文件，实际上写入的是线单元
	/// @param _name 文件名称
	void writeGeometryVtkFile(std::string _name);

	/// @brief 设置重分网区域
	/// @param _deformedlist 重分网区域列表
	virtual void setDeformedDomain(const std::vector<int> _deformedlist);

	/// @brief 设置运动区域
	/// @param _movingmap 运动区域到运动属性的映射关系
	virtual void setMovingPart(const std::map<int, FEMMovingPart*> _movingmap);
	void setMeshManager(FEMMeshManager* _meshmanager);

	/// @brief 更新电流（动态特性计算用到）
	/// @param domain 几何区域
	/// @param current 电流大小
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
	std::vector<int> deformedlist;	//形变区域
	std::map<int, FEMMovingPart*>movingmap;	//运动区域

	int maxitersteps;
	double maxerror;

	int num_freenodes;	//自由节点数目
	std::vector<int> node_reorder;	//前num_dof个元素对应非边界节点，之后的元素对应第一类边界条件
	std::vector<int> node_pos;	//原节点编号对应的reorder后的节点编号

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

