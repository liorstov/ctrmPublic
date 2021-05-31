#pragma once
#include <vector>

class Geophone
{
public:
	Geophone(int, float, float, float);
	int index;
	float x, y, z;
	std::vector<float> U;
};

