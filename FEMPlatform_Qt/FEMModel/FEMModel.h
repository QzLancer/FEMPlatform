#pragma once

#include <vector>
#include <map>
#include <string>

#include "../FEMMaterial.h"
#include "../FEMBoundary.h"
#include "../FEMMovingPart.h"

#define PI 3.14159265358979323846

//模型类，需要使用新模型时，继承该类/创建该类对象，然后设置好private中的全部必要参数
//相比于类更像是一个结构体
class FEMModel 
{
protected:
	FEMModel();
	/// <summary>
	/// 设置模型名称
	/// </summary>
	virtual void setModelName() = 0;

	/// <summary>
	/// 设置维度
	/// </summary>
	virtual void setDimension() = 0;

	/// <summary>
	/// 设置几何or网格文件名称
	/// </summary>
	virtual void setFile() = 0;	

	/// <summary>
	/// 添加非线性电磁材料（BH曲线）
	/// </summary>
	/// <param name="_name">材料名称</param>
	/// <param name="_bhpoints">BH曲线插值点的数量</param>
	/// <param name="_bdata">B数组</param>
	/// <param name="_hdata">H数组</param>
	void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata);

	/// <summary>
	/// 添加线性电磁材料
	/// </summary>
	/// <param name="_name">材料名称</param>
	/// <param name="_mu">相对磁导率</param>
	void addLinearMaterial(std::string _name, double _mu);

	/// @brief 添加永磁材料
	/// @param _name 材料名称
	/// @param _mu 相对磁导率
	/// @param _h_c 矫顽力
	/// @param _b_r 剩磁
	/// @param _theta_m 永磁相对x轴正方向的角度
	void addPermMagentMaterial(std::string _name, double _mu, double _h_c, double _b_r, double _theta_m);

	/// @brief 添加线圈
	/// @param _name 材料名称
	/// @param _coil 线圈类
	void addCoil(std::string _name, FEMCoil _coil);

	/// @brief 设置单元材料类型
	virtual void createElement2Material() = 0;

	/// @brief 设置负载
	virtual void bulidGeometry2Load() = 0;

	/// @brief 设置边界条件
	virtual void buildGeometry2Constrain() = 0;

	/// @brief 设置比例系数（Gmsh网格和COMSOL网格的尺寸才是不同）
	virtual void setUnitRatio() = 0; 

	/// @brief 设置形变区域
	virtual void buildGeometry2Deformed() = 0;

	/// @brief 设置运动区域
	virtual void buildGeometry2MovingPart() = 0;

public:
	/// <summary>
	/// 对模型进行初始化，被调用的函数在主函数中为虚函数，此处使用Template模式
	/// </summary>
	void init();

	enum class DIMENSION {
		ONE = 1,
		D2AXISM = 2,
		THREE = 3,
		D2PLANE = 4
	};

	/// <summary>
	/// 获得求解问题的维度
	/// </summary>
	/// <returns>1：一维；2：二维轴对称；3：三维；4：二维平面</returns>
	DIMENSION getDimension() const;

	/// <summary>
	/// </summary>
	/// <returns>COMSOL网格文件的完整路径及名称</returns>
	std::string getMeshFile() const;

	/// <summary>
	/// </summary>
	/// <returns>用于Gmsh分网的.geo几何文件名称</returns>
	std::string getGeoFile() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>材料列表</returns>
	std::vector<FEMMaterial*> getMaterialList() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>几何域到材料之间的映射关系</returns>
	std::map<int, FEMMaterial*> getMaterialMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>激励区域的电流大小</returns>
	std::map<int, double> getLoadMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>边界区域到边界条件类型的映射关系</returns>
	std::map<int, FEMBoundary*> getBoundaryMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>模型名称</returns>
	std::string getModelName() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>网格尺寸比例系数（Gmsh和COMSOL的网格尺寸使用的是不同的单位，需要换算）</returns>
	double getUnitRatio() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>重分网区域的编号</returns>
	std::vector<int> getDeformedList() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>运动区域的编号</returns>
	std::map<int, FEMMovingPart*> getMovingMap() const;

protected:
	/// @brief 默认单位是m，如果输入的几何参数是mm，则需要将ratio置为0.001；
	double unitratio{1};
	/// @brief 模型名称
	std::string modelname;
	/// @brief 维度
	DIMENSION dimension;
	/// @brief 几何文件名称
	std::string geofile;
	/// @brief 网格文件名称
	std::string meshfile;
	/// @brief 材料列表
	std::vector<FEMMaterial*> materiallist;
	/// @brief 材料到几何区域的映射关系
	std::map<int, FEMMaterial*> materialmap;
	/// @brief 电流大小到几何区域的映射关系
	std::map<int, double> loadmap;
	/// @brief 边界条件到几何区域的映射关系
	std::map<int, FEMBoundary*> boundarymap;
	/// @brief 形变区域列表
	std::vector<int> deformedlist;
	/// @brief 运动条件到运动区域的映射关系
	std::map<int, FEMMovingPart*>movingmap;
};