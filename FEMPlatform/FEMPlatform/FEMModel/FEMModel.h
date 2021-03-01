#pragma once

#include <vector>
#include <map>
#include <string>

#include "../FEMMaterial.h"
#include "../FEMBoundary.h"

#define PI 3.14159265358979323846

//ģ���࣬��Ҫʹ����ģ��ʱ���̳и���/�����������Ȼ�����ú�private�е�ȫ����Ҫ����
//������������һ���ṹ��
class FEMModel 
{
protected:
	FEMModel();
	virtual void setModelName() = 0;
	virtual void setdimension() = 0;
	virtual void setMeshFile() = 0;	/*���������ļ�����*/
	void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata);/*��ӷ����Բ���*/
	void addLinearMaterial(std::string _name, double _mu);/*������Բ���*/
	virtual void createElement2Material() = 0;	/*���õ�Ԫ��������*/
	virtual void bulidGeometry2Load() = 0;		/*���ø���*/
	virtual void buildGeometry2Constrain() = 0;	/*���ñ߽�����*/

public:
	void init();	/*���ģʽ��Templateģʽ*/
	enum class DIMENSION {
		ONE = 1,
		TWO = 2,
		THREE = 3
	};
	DIMENSION getDimension() const;
	std::string getMeshFile() const;
	std::vector<FEMMaterial*> getMaterialList() const;
	std::map<int, FEMMaterial*> getMaterialMap() const;
	std::map<int, double> getLoadMap() const;
	std::map<int, FEMBoundary*> getBoundaryMap() const;
	std::string getModelName() const;


protected:
	std::string modelname;
	DIMENSION dimension;
	std::string meshfile;
	std::vector<FEMMaterial*> materiallist;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
};