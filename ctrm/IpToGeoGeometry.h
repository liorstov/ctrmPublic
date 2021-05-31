#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
struct DeltaVelocity
{
	DeltaVelocity(float, float);
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
	IpToGeoGeometry(int _index, float _surface, float _betta, float _verD);
	void calculateTimeDelta(std::vector<std::pair<float, float>> const& velocities);
	float calcAvarageVelo(std::vector<std::pair<float, float>> const& velocities, float begin, float end);
	std::string writeInfo(int IpIndex);
	float surfaceDist, bettaAngle, verDepth;
	int GeoIndex;
	std::vector<RadiusTimeDelta> RadiusTimeValues;
	float dt = 0.500000;

};



