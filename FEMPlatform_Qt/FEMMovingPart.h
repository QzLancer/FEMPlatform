#pragma once
struct MovingLimit
{
	double min{ 0 };
	double max{ 0 };
};

/// @brief 运动区域类，
class FEMMovingPart
{
public:
	FEMMovingPart();
	~FEMMovingPart();
	
	/// @brief 运动方向，通过长度为3的数组表示，三个元素分别表示x, y, z方向，将对应方向的值置为1，表示组件运动方向
	int direction[3]{ 0, 0, 0 };

	/// @brief 位移限制，每个方向(xyz)的运动范围限制由min和max决定
	MovingLimit limit[3];

	/// @brief 
	/// @return 运动组件的质量
	double getMass() const;

	/// @brief 设置运动组件的质量
	/// @param _mass 运动组件的质量
	void setMass(const double _mass);

	/// @brief 
	/// @return 用于弹簧力计算的数组列表的元素数目
	int getForceNodeSize() const;

	/// @brief 设置弹簧力数组列表（位置-弹簧力大小的映射关系）
	/// @param _nodesize 插值节点数目
	/// @param _pos 位置数组
	/// @param _force 弹簧力数组
	void setSpringForce(int _nodesize, double* _pos, double* _force);

	/// @brief 通过插值的方式计算弹簧力
	/// @param _pos 运动部件所在的位置
	/// @return 弹簧力大小
	double getSpringForce(const double _pos) const;

private:
	/// @brief 弹簧力数组列表元素数目
	int forcenodesize;
	/// @brief 运动部件质量
	double mass;
	/// @brief 运动部件位置列表
	double* position;
	/// @brief 弹簧力列表
	double* force;
};
