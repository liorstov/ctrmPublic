#pragma once
#include <vector>
#include <tuple>
#include <algorithm>
#include <stdio.h>

// velocity for each radius
struct velocity {
	velocity(int);
	int value;
	float semb;

};
// each sample has a vector of multipass velocities
struct sampleData {
	sampleData(int, int ,int);
	int sampN;
	float semblance;
	std::vector<velocity> velo;
};

//each imagePOint has a vector of samples 
struct ImageP
{
	ImageP(int, float, float, float, int ,int,int  ,int,int); 
	float x, y, z;
	int IPindex;
	std::vector<std::tuple<int,float,float>> timeDeltas;
	std::vector<sampleData> samples;

};