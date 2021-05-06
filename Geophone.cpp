#include "Geophone.h"

Geophone::Geophone(int _ind,int _x, int _y, int _z)
{
	index = _ind;
	x = _x;
	y = _y;
	z = _z;

	U.reserve(sizeof(float) * 500);
}
