#pragma once
struct MovingLimit
{
	double min{ 0 };
	double max{ 0 };
};

class FEMMovingPart
{
public:
	//�˶�����ͨ������Ϊ3�������ʾ������Ԫ�طֱ��ʾx, y, z����
	int direction[3]{ 0, 0, 0 };

	//λ�����ƣ�ͨ������Ϊ6�������ʾ��ÿ��������min��max����
	MovingLimit limit[3];

};
