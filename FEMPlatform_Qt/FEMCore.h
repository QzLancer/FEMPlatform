#pragma once
#include "MeshManager/FEMMeshManager.h"
#include "MeshManager/FEM2DMeshManager.h"
#include "MeshManager/FEM3DMeshManager.h"

#include "FEMSolver/FEMSolver.h"
#include "FEMSolver/FEM2DSolver.h"
#include "FEMSolver/FEM2DNRSolver.h"
#include "FEMSolver/FEM2DNDDRSolver.h"
#include "FEMSolver/FEM2DNDDRCUDASolver.cuh"
#include "FEMSolver/FEM2DSchwarzSolver.h"

#include "MatrixSolver/SluMTMatrixSolver.h"

#include "FEMModel/FEMModel.h"

#include <iostream>
#include <string>
using namespace std;

/// @brief 提供用户操作的接口，给出用户在图形界面可以进行的设置
class FEMCore
{
public:
	FEMCore();
	~FEMCore();

	/// @brief 将模型中的几何文件名称或网格文件名称传递给FEMMeshManager
	/// @param _model 
	void setModel(FEMModel* _model);

	/// @brief 选择求解静态特性还是动态特性，目前该函数为空，静动态特性在solve()中进行判断
	/// @param analysistype 分析类型，static or dynamic
	void setAnalysisType(string analysistype);

	/// @brief 配置求解算法
	/// @param strategy NR, NDDR, NDDRGPU or DD
	void setFEMSolveStrategy(string strategy);

	/// @brief 设置矩阵求解器
	/// @param matrixsolver SuperLU_MT
	void setMatrixSolver(string matrixsolver);

	/// @brief FEM求解，静态特性或者动态特性
	/// @param analysistype 分析类型，static or dynamic
	void solve(string analysistype);

	/// @brief 后处理，输出vtk文件，通过paraview后处理
	void postprocess();

	/// @brief 设置最大迭代步数（需要改进，哪一层迭代的最大步数）
	/// @param _maxitersteps 
	void setMaxIterSteps(const int _maxitersteps);

	/// @brief 设置迭代允许的最大误差（需要改进，哪一层迭代的最大误差）
	/// @param _error 
	void setMaxError(const double _error);

private:
	FEMModel* model;
	FEMMeshManager* meshmanager;
	FEMSolver* solver;
};

