#pragma once
#include <vector>
#include <tuple>
#include <algorithm>
#include "IpToGeoGeometry.h"
#include <stdio.h>

class velocity {
public:
	velocity(int);
	int value;
	float semb;

};
struct sampleData {
	sampleData(int, int ,int);
	int sampN;
	float semblance;
	std::vector<velocity> velo;
};


class ImageP
{
public:
	ImageP(int, float, float, float, int ,int,int  ,int,int); 
	float x, y, z;
	int IPindex;
	std::vector<std::tuple<int,float,float>> timeDeltas;
	std::vector<sampleData> samples;

};