#pragma once
#include <vector>
#include <tuple>
#include <algorithm>
#include "IpToGeoGeometry.h"
class ImageP
{
public:
	ImageP(int,float, float, float);
	float x, y, z;
	int IPindex;
	std::vector<IpToGeoGeometry> IpGeoCalc;
	std::vector<std::tuple<int,int,float>> SampleSemblance;


};

