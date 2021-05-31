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
#include <Eigen/Dense>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace std;
class box
{
	std::vector<Geophone> vGeophones;
	std::vector<ImageP> vImagePoints;

public:
	box(int, int,int,int,int,int,int,int, int, int);
	void readCoord(string file);
	void setCoord(Eigen::MatrixXf);
	void setEnergy(Eigen::MatrixXf);
	void readVelo(string file);
	void readEnergy(string file);
	Eigen::MatrixXd getEnergy();
	void createImageSpace();
	int getIP(int point, int sample);
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
	int xmax = 5;
	int ymax = 5;
	int dxTrace = 1;
	int dyLines = 10;
	int coordcy0 = 45;
	int nsamp = 1000;
		
	int jbeg = 0;
	int jend = 700;
	int startRadius = 1;
	int endRadius = 5;
	float rmin = 1;
	float rmax = 35;
	float dr = 1;
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
	
};

