#include "Geophone.h"

Geophone::Geophone(int _ind, float _x, float _y, float _z)
{
	index = _ind;
	x = _x;
	y = _y;
	z = _z;

	U.reserve(sizeof(float) * 500);
}
