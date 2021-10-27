#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
struct DeltaVelocity
{
	DeltaVelocity(float, float,float);
	float referenceVelocity {0};
	float actualVelocity{ 0 };
	float timeDelta{ 0 };
};
struct RadiusTimeDelta
{
	RadiusTimeDelta(float _radius);
	
	float radius;
	float distanceToIP;
	std::vector<DeltaVelocity> VelocityRangeTimeDelta;
};
class IpToGeoGeometry
{
public:
	IpToGeoGeometry(int _index, float _horizontal,float _vertival, float _distance, float xdist,float ydist);
	void calculateTimeDelta(std::vector<std::pair<float, float>> const& velocities);
	float calcAvarageVelo(std::vector<std::pair<float, float>> const& velocities, float begin, float end);
	std::string writeInfo(int IpIndex);
	void RP();
	float horizontal, vertical, distance,xdist,ydist ,RadPatternPwave , RadPatternSwave;
	int GeoIndex;
	std::vector<DeltaVelocity> VelocityRangeTimeDelta;
	float dt = 0.500000;

};



