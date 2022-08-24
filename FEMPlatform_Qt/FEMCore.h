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

/// @brief �ṩ�û������Ľӿڣ������û���ͼ�ν�����Խ��е�����
class FEMCore
{
public:
	FEMCore();
	~FEMCore();

	/// @brief ��ģ���еļ����ļ����ƻ������ļ����ƴ��ݸ�FEMMeshManager
	/// @param _model 
	void setModel(FEMModel* _model);

	/// @brief ѡ����⾲̬���Ի��Ƕ�̬���ԣ�Ŀǰ�ú���Ϊ�գ�����̬������solve()�н����ж�
	/// @param analysistype �������ͣ�static or dynamic
	void setAnalysisType(string analysistype);

	/// @brief ��������㷨
	/// @param strategy NR, NDDR, NDDRGPU or DD
	void setFEMSolveStrategy(string strategy);

	/// @brief ���þ��������
	/// @param matrixsolver SuperLU_MT
	void setMatrixSolver(string matrixsolver);

	/// @brief FEM��⣬��̬���Ի��߶�̬����
	/// @param analysistype �������ͣ�static or dynamic
	void solve(string analysistype);

	/// @brief �������vtk�ļ���ͨ��paraview����
	void postprocess();

	/// @brief ������������������Ҫ�Ľ�����һ��������������
	/// @param _maxitersteps 
	void setMaxIterSteps(const int _maxitersteps);

	/// @brief ���õ���������������Ҫ�Ľ�����һ������������
	/// @param _error 
	void setMaxError(const double _error);

private:
	FEMModel* model;
	FEMMeshManager* meshmanager;
	FEMSolver* solver;
};

