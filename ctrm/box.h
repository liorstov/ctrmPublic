#pragma once
#include <vector>
#include <tuple>
#include <iostream>
#include <string>
#include <stdio.h>
#include "Geophone.h"
#include "ImageP.h"
#include "IpToGeoGeometry.h"
#include <fstream>
#include <cmath>
#include <chrono>
#include "progressbar.hpp"
#include <algorithm>
#include <thread>
#include <omp.h>
#include "npy.hpp"


using namespace std;
class box
{
	std::vector<Geophone> vGeophones;
	std::vector<ImageP> vImagePoints;

public:
	box(int traces, int lines, int tracesMin, int linesMin, int startRad, int endRad, int _dx, int _dy, int _jbeg, int _jend, int _vrange, int _dv, int _minDist, int _windowSize, int _dr);
	void readCoord(string file);
	
	void readEnergy(string file);
	void readVelo(string file);
	void createImageSpace();
	std::vector<int> getCoord();
	float getGeophoneEnergy(int , int );

	float getGeoZ(int z, int y);

	void CalcSurfaceDist();
	void CalcTimeDeltaOnly();

	//void corrolationOnGeo(bool RP);
	float calcAvarageVelo(float , float);
	//Eigen::MatrixXd  getIP();
	//void writeIP();
	void writeSemblence();
	void writeSemblenceNpy(std::string file);
	void writeTimeDeltasNpy(std::string file);
	std::vector<pair<float,float>> vRadiusVelo;
	float highestEnergy{ 0.0f };
	float minimumDepth = 1000000;
	float ntrace = 72;
	float numrecl = 1;
	float numsh = 1;
	float numshli = 1;
	float xmax = 5;
	float xmin = 5;
	float ymin = 5;
	float ymax = 5;
	int dxTrace = 1;
	int dyLines = 10;
	int nsamp=1000;
	int windowSize;
	int jbeg = 400;
	int jend = 700;
	int minDist = 1000;
	int startRadius = 1;
	int endRadius = 5;
	float rmin = 1;
	float rmax = 35;
	int dr = 1;
	float v0 = 500;
	float dt = (0.500000f)/1000.0f;
	float resamp = 1;
	float vmin = 9999;
	float vmax = 9999;
	int dv = 30;
	float itim1 = 9999;
	float itim2 = 9999;
	float offmin = 0;
	float offmax = 999;
	int vRange = 150;
	int signalPosition = 0;
	
};

