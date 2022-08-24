#pragma once

#include <string>

/// @brief �߽������࣬Ŀǰֻ֧����ԳƱ߽�����ɱ߽�����������Ԫ�����в��ý��ڵ����ɶ�ֱ�Ӹ�ֵΪ0�ļ��㷽����
class FEMBoundary
{
public:
	virtual std::string getName() const = 0;
	int getBoundaryType() const;

protected:
	int boundarytype;	//�߽�����������ȷ�����ɶ�

};

class FEMAxialSymmetry :
	public FEMBoundary 
{
public:
	FEMAxialSymmetry();
	virtual std::string getName() const override;
};

class FEMMagneticInsulation :
	public FEMBoundary
{
public:
	FEMMagneticInsulation();
	virtual std::string getName() const override;
};