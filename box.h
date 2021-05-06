#pragma once
#include <vector>
#include <tuple>
#include <iostream>
#include <string>
#include "Geophone.h"
#include "ImageP.h"
#include "IpToGeoGeometry.h"
#include <fstream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
using namespace std;
class box
{
	std::vector<Geophone> vGeophones;
	std::vector<ImageP> vImagePoints;

public:
	box(int, int,int,int);
	void readCoord(string file);
	void readVelo(string file);
	void readEnergy(string file);
	void createImageSpace();
	float getIP(int point, int sample);
	std::vector<int> getCoord();
	float getGeophoneEnergy(int , int );

	float getGeoZ(int z, int y);

	void CalcSurfaceDist();

	void corrolationOnGeo();
	float calcAvarageVelo(float begin, float end);
	Eigen::MatrixXd  getSample();
	void writeIP();
	void writeSemblence();
	std::vector<pair<float,float>> vRadiusVelo;
	
	float minimumDepth = FLT_MAX;
	float ntrace = 72;
	float numrecl = 1;
	float numsh = 1;
	float numshli = 1;
	int numcen = 5;
	int numlin = 5;
	float dxc = 1.000000;
	float dyc = 10.000000;
	float coordcy0 = 45.000000;
	int nsamp = 500;
		
	int jbeg = 0;
	int jend = 400;
	int startRadius = 1;
	int endRadius = 5;
	float rmin = 1;
	float rmax = 35;
	float dr = 1;
	float v0 = 500;
	float xmin = 0;
	float xmax = 999;
	float dt = (0.500000f)/1000.0f;
	float resamp = 1;
	float vmin = 9999;
	float vmax = 9999;
	float dv = 20;
	float itim1 = 9999;
	float itim2 = 9999;
	float offmin = 0;
	float offmax = 999;
	float vRange = 100;
	
};

