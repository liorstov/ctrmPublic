#include "IpToGeoGeometry.h"

IpToGeoGeometry::IpToGeoGeometry(int _index,float _surface, float _betta, float _verD)
{
	GeoIndex = _index;
	surfaceDist = _surface;
	bettaAngle = _betta;
	verDepth = _verD;
	RadiusTimeValues.reserve(sizeof(RadiusTimeDelta) * 100);
}

void IpToGeoGeometry::calculateTimeDelta(std::vector<std::pair<float,float>> const& velocities)
{

	float VV,radius,CurrectVelocity, deltaTime;

	for (auto const& value : velocities) {
		radius = value.first;
		CurrectVelocity = value.second;
		VV = calcAvarageVelo(velocities, verDepth, radius);
		deltaTime = sqrtf(powf(surfaceDist, 2) + powf(verDepth + radius, 2) - 2 * sinf(bettaAngle) * surfaceDist * verDepth / VV - verDepth / CurrectVelocity);
		//VeloTimeDelta.push_back(deltaTime / dt);
	}
}

float IpToGeoGeometry::calcAvarageVelo(std::vector<std::pair<float, float>> const& velocities, float begin, float end)
{
	float radius, CurrectVelocity, avarage = 0;
	for (auto const& value : velocities) 
	{
		radius = value.first;
		CurrectVelocity = value.second;

		if (radius >= begin && radius <= end) {
			avarage += (1 / CurrectVelocity);
		}
	}

	return((end - begin) / avarage);
}

std::string IpToGeoGeometry::writeInfo(int IpIndex)
{
	char res[100];
	sprintf_s(res,100, "%d, %1.0f, %1.0f, %1.0f",  GeoIndex,surfaceDist,bettaAngle,verDepth);
	return(res);
}

RadiusTimeDelta::RadiusTimeDelta(float _radius) : radius(_radius)
{
}

DeltaVelocity::DeltaVelocity(float _actualVelo, int _timeDelta) : actualVelocity(_actualVelo), timeDelta(_timeDelta)
{
}
