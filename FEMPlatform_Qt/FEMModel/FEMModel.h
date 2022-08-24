#pragma once

#include <vector>
#include <map>
#include <string>

#include "../FEMMaterial.h"
#include "../FEMBoundary.h"
#include "../FEMMovingPart.h"

#define PI 3.14159265358979323846

//ģ���࣬��Ҫʹ����ģ��ʱ���̳и���/�����������Ȼ�����ú�private�е�ȫ����Ҫ����
//������������һ���ṹ��
class FEMModel 
{
protected:
	FEMModel();
	/// <summary>
	/// ����ģ������
	/// </summary>
	virtual void setModelName() = 0;

	/// <summary>
	/// ����ά��
	/// </summary>
	virtual void setDimension() = 0;

	/// <summary>
	/// ���ü���or�����ļ�����
	/// </summary>
	virtual void setFile() = 0;	

	/// <summary>
	/// ��ӷ����Ե�Ų��ϣ�BH���ߣ�
	/// </summary>
	/// <param name="_name">��������</param>
	/// <param name="_bhpoints">BH���߲�ֵ�������</param>
	/// <param name="_bdata">B����</param>
	/// <param name="_hdata">H����</param>
	void addNonlinearMaterial(std::string _name, int _bhpoints, double* _bdata, double* _hdata);

	/// <summary>
	/// ������Ե�Ų���
	/// </summary>
	/// <param name="_name">��������</param>
	/// <param name="_mu">��Դŵ���</param>
	void addLinearMaterial(std::string _name, double _mu);

	/// @brief ������Ų���
	/// @param _name ��������
	/// @param _mu ��Դŵ���
	/// @param _h_c ������
	/// @param _b_r ʣ��
	/// @param _theta_m �������x��������ĽǶ�
	void addPermMagentMaterial(std::string _name, double _mu, double _h_c, double _b_r, double _theta_m);

	/// @brief �����Ȧ
	/// @param _name ��������
	/// @param _coil ��Ȧ��
	void addCoil(std::string _name, FEMCoil _coil);

	/// @brief ���õ�Ԫ��������
	virtual void createElement2Material() = 0;

	/// @brief ���ø���
	virtual void bulidGeometry2Load() = 0;

	/// @brief ���ñ߽�����
	virtual void buildGeometry2Constrain() = 0;

	/// @brief ���ñ���ϵ����Gmsh�����COMSOL����ĳߴ���ǲ�ͬ��
	virtual void setUnitRatio() = 0; 

	/// @brief �����α�����
	virtual void buildGeometry2Deformed() = 0;

	/// @brief �����˶�����
	virtual void buildGeometry2MovingPart() = 0;

public:
	/// <summary>
	/// ��ģ�ͽ��г�ʼ���������õĺ�������������Ϊ�麯�����˴�ʹ��Templateģʽ
	/// </summary>
	void init();

	enum class DIMENSION {
		ONE = 1,
		D2AXISM = 2,
		THREE = 3,
		D2PLANE = 4
	};

	/// <summary>
	/// �����������ά��
	/// </summary>
	/// <returns>1��һά��2����ά��Գƣ�3����ά��4����άƽ��</returns>
	DIMENSION getDimension() const;

	/// <summary>
	/// </summary>
	/// <returns>COMSOL�����ļ�������·��������</returns>
	std::string getMeshFile() const;

	/// <summary>
	/// </summary>
	/// <returns>����Gmsh������.geo�����ļ�����</returns>
	std::string getGeoFile() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>�����б�</returns>
	std::vector<FEMMaterial*> getMaterialList() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>�����򵽲���֮���ӳ���ϵ</returns>
	std::map<int, FEMMaterial*> getMaterialMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>��������ĵ�����С</returns>
	std::map<int, double> getLoadMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>�߽����򵽱߽��������͵�ӳ���ϵ</returns>
	std::map<int, FEMBoundary*> getBoundaryMap() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>ģ������</returns>
	std::string getModelName() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>����ߴ����ϵ����Gmsh��COMSOL������ߴ�ʹ�õ��ǲ�ͬ�ĵ�λ����Ҫ���㣩</returns>
	double getUnitRatio() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>�ط�������ı��</returns>
	std::vector<int> getDeformedList() const;

	/// <summary>
	/// 
	/// </summary>
	/// <returns>�˶�����ı��</returns>
	std::map<int, FEMMovingPart*> getMovingMap() const;

protected:
	/// @brief Ĭ�ϵ�λ��m���������ļ��β�����mm������Ҫ��ratio��Ϊ0.001��
	double unitratio{1};
	/// @brief ģ������
	std::string modelname;
	/// @brief ά��
	DIMENSION dimension;
	/// @brief �����ļ�����
	std::string geofile;
	/// @brief �����ļ�����
	std::string meshfile;
	/// @brief �����б�
	std::vector<FEMMaterial*> materiallist;
	/// @brief ���ϵ����������ӳ���ϵ
	std::map<int, FEMMaterial*> materialmap;
	/// @brief ������С�����������ӳ���ϵ
	std::map<int, double> loadmap;
	/// @brief �߽����������������ӳ���ϵ
	std::map<int, FEMBoundary*> boundarymap;
	/// @brief �α������б�
	std::vector<int> deformedlist;
	/// @brief �˶��������˶������ӳ���ϵ
	std::map<int, FEMMovingPart*>movingmap;
};