#pragma once
#include <vector>
#include <tuple>
#include <iostream>
#include <string>
#include <stdio.h>
#include "Geophone.h"
#include "ImageP.h"
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
	box(int traces, int lines, int tracesMin, int linesMin, int startRad, int endRad, int _dx, int _dy, int _jbeg, int _jend, int _vrange, int _dv, int _minDist, int _windowSize, int _dr,float _numOfThreadsPercent, float _sampleRate);
	void readCoord(string file);
	void readCoordnumpy(string file);
	
	void readEnergy(string file);
	void readEnergyFromNpy(string file);
	void readVelo(string file);
	void createImageSpace();
	std::vector<int> getCoord();

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
	float highestEnergy;
	float minimumDepth;
	float ntrace;
	float numrecl;
	float numsh;
	float numshli;
	float xmax;
	float xmin;
	float ymin;
	float ymax;
	int dxTrace;
	int dyLines;
	int nsamp;
	int windowSize;
	int jbeg;
	int jend;
	int minDist;
	int startRadius;
	int endRadius;
	float rmin;
	float rmax;
	int dr;
	float v0;
	float dt;
	float resamp;
	float vmin;
	float vmax;
	int dv;
	float itim1;
	float itim2;
	float offmin;
	float offmax;
	int vRange;
	int signalPosition;
	float numOfThreadsPercent;
};

