#include "FEMMovingPart.h"

FEMMovingPart::FEMMovingPart() :
	forcenodesize(0),
	position(nullptr),
	force(nullptr)
{
}

FEMMovingPart::~FEMMovingPart()
{
	if (position != nullptr) {
		delete[] position;
	}
	if (force != nullptr) {
		delete[] force;
	}
}

double FEMMovingPart::getMass() const
{
	return mass;
}

void FEMMovingPart::setMass(const double _mass)
{
	mass = _mass;
}

int FEMMovingPart::getForceNodeSize() const
{
	return forcenodesize;
}

void FEMMovingPart::setSpringForce(int _nodesize, double* _pos, double* _force)
{
	forcenodesize = _nodesize;
	position = _pos;
	force = _force;
}

double FEMMovingPart::getSpringForce(const double _pos) const
{
	if (_pos > position[forcenodesize - 1]) {
		int  i = forcenodesize - 1;
		double k = (force[i] - force[i - 1]) / (position[i] - position[i - 1]);
		return force[i] + k * (_pos - position[i]);
	}
	else if (_pos < position[0]) {
		return force[0];
	}
	else
	{
		for (int i = 0; i < forcenodesize - 1; i++) {
			if (_pos >= position[i] && _pos <= position[i + 1]) {
				double k = (force[i + 1] - force[i]) / (position[i + 1] - position[i]);
				return force[i] + k * ((_pos - position[i]));
			}
		}
	}
}
