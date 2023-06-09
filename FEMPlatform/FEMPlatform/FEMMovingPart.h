#pragma once
struct MovingLimit
{
	double min{ 0 };
	double max{ 0 };
};

class FEMMovingPart
{
public:
	//运动方向，通过长度为3的数组表示，三个元素分别表示x, y, z方向
	int direction[3]{ 0, 0, 0 };

	//位移限制，通过长度为6的数组表示，每个方向由min和max决定
	MovingLimit limit[3];

};
