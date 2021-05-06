#pragma once
#include <vector>

class Geophone
{
public:
	Geophone(int,int, int, int);
	int index, x, y, z;
	std::vector<float> U;
};

