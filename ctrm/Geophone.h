#pragma once
#include <vector>


// Geophone class
class Geophone
{
public:
	Geophone(int, float, float, float);
	int index;
	float x, y, z;
	std::vector<float> U;
};

