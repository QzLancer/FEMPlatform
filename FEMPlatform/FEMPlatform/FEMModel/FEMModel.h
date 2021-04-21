#pragma once

#include <vector>
#include <map>
#include <string>

#include "../FEMMaterial.h"
#include "../FEMBoundary.h"

#define PI 3.14159265358979323846

//模型类，需要使用新模型时，继承该类/创建该类对象，然后设置好private中的全部必要参数
//相比于类更像是一个结构体
class FEMModel 
{
protected:
	FEMModel();
	virtual void setModelName() = 0;
	virtual void setDimension() = 0;
	virtual void setFile() = 0;	/*设置几何or网格文件名称*/
	void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata);	/*添加非线性材料*/
	void addLinearMaterial(std::string _name, double _mu);	  /*添加线性材料*/
	void addPermMagentMaterial(std::string _name, double _mu, double _h_c, double _b_r, double _theta_m);
	void addCoil(std::string _name, FEMCoil _coil);
	virtual void createElement2Material() = 0;	/*设置单元材料类型*/
	virtual void bulidGeometry2Load() = 0;		/*设置负载*/
	virtual void buildGeometry2Constrain() = 0;	/*设置边界条件*/
	virtual void setUnitRatio() = 0; 

public:
	void init();	/*设计模式：Template模式*/
	enum class DIMENSION {
		ONE = 1,
		TWO = 2,
		THREE = 3
	};
	DIMENSION getDimension() const;
	std::string getMeshFile() const;
	std::string getGeoFile() const;
	std::vector<FEMMaterial*> getMaterialList() const;
	std::map<int, FEMMaterial*> getMaterialMap() const;
	std::map<int, double> getLoadMap() const;
	std::map<int, FEMBoundary*> getBoundaryMap() const;
	std::string getModelName() const;
	double getUnitRatio() const;


protected:
	double unitratio{1};	//默认单位是m，如果输入的几何参数是mm，则需要将ratio置为0.001；
	std::string modelname;
	DIMENSION dimension;
	std::string geofile;
	std::string meshfile;
	std::vector<FEMMaterial*> materiallist;
	std::map<int, FEMMaterial*> materialmap;
	std::map<int, double> loadmap;
	std::map<int, FEMBoundary*> boundarymap;
};