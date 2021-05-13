#pragma once
struct MovingLimit
{
	double min{ 0 };
	double max{ 0 };
};

class FEMMovingPart
{
public:
	FEMMovingPart();
	~FEMMovingPart();
	//�˶�����ͨ������Ϊ3�������ʾ������Ԫ�طֱ��ʾx, y, z����
	int direction[3]{ 0, 0, 0 };

	//λ�����ƣ�ͨ������Ϊ6�������ʾ��ÿ��������min��max����
	MovingLimit limit[3];

	double getMass() const;
	void setMass(const double _mass);
	int getForceNodeSize() const;
	void setSpringForce(int _nodesize, double* _pos, double* _force);
	double getSpringForce(const double _pos) const;

private:
	int forcenodesize;
	double mass;
	double* position;
	double* force;
};
