#pragma once
struct MovingLimit
{
	double min{ 0 };
	double max{ 0 };
};

/// @brief �˶������࣬
class FEMMovingPart
{
public:
	FEMMovingPart();
	~FEMMovingPart();
	
	/// @brief �˶�����ͨ������Ϊ3�������ʾ������Ԫ�طֱ��ʾx, y, z���򣬽���Ӧ�����ֵ��Ϊ1����ʾ����˶�����
	int direction[3]{ 0, 0, 0 };

	/// @brief λ�����ƣ�ÿ������(xyz)���˶���Χ������min��max����
	MovingLimit limit[3];

	/// @brief 
	/// @return �˶����������
	double getMass() const;

	/// @brief �����˶����������
	/// @param _mass �˶����������
	void setMass(const double _mass);

	/// @brief 
	/// @return ���ڵ���������������б��Ԫ����Ŀ
	int getForceNodeSize() const;

	/// @brief ���õ����������б�λ��-��������С��ӳ���ϵ��
	/// @param _nodesize ��ֵ�ڵ���Ŀ
	/// @param _pos λ������
	/// @param _force ����������
	void setSpringForce(int _nodesize, double* _pos, double* _force);

	/// @brief ͨ����ֵ�ķ�ʽ���㵯����
	/// @param _pos �˶��������ڵ�λ��
	/// @return ��������С
	double getSpringForce(const double _pos) const;

private:
	/// @brief �����������б�Ԫ����Ŀ
	int forcenodesize;
	/// @brief �˶���������
	double mass;
	/// @brief �˶�����λ���б�
	double* position;
	/// @brief �������б�
	double* force;
};
