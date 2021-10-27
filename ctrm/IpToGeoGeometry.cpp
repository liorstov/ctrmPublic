#include "IpToGeoGeometry.h"

IpToGeoGeometry::IpToGeoGeometry(int _index, float _horizontal, float _vertical, float _distance, float xdist, float ydist):GeoIndex(_index), horizontal(_horizontal), vertical(_vertical), distance(_distance), xdist(xdist),ydist(ydist)
{
	this->RP();
}

void IpToGeoGeometry::calculateTimeDelta(std::vector<std::pair<float,float>> const& velocities)
{  
	
	//float VV,radius,CurrectVelocity, deltaTime;

	//for (auto const& value : velocities) {
	//	radius = value.first;
	//	CurrectVelocity = value.second;
	//	VV = calcAvarageVelo(velocities, verDepth, radius);
	//	deltaTime = sqrtf(powf(surfaceDist, 2) + powf(verDepth + radius, 2) - 2 * sinf(bettaAngle) * surfaceDist * verDepth / VV - verDepth / CurrectVelocity);
	//	//VeloTimeDelta.push_back(deltaTime / dt);
	//}
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
	sprintf_s(res,100, "%d, %1.0f, %1.0f, %1.0f",  GeoIndex,this->horizontal, this->vertical,this->distance);
	return(res);
}

void IpToGeoGeometry::RP()
{
	float cosPhi, sinPhi, cosTheta, sinTheta;
	cosPhi = -this->ydist / this->horizontal;
	sinPhi = -this->xdist / this->horizontal;
	cosTheta = -this->ydist / sqrtf(powf(this->ydist,2)+powf(this->vertical,2));
	sinTheta = -this->vertical / sqrtf(powf(this->ydist,2)+powf(this->vertical,2));
	sinTheta *= -this->ydist / fabsf(this->ydist);

	this->RadPatternPwave = -sinTheta * cosTheta * cosPhi/this->distance;
	this->RadPatternSwave = -this->RadPatternPwave;
}

RadiusTimeDelta::RadiusTimeDelta(float _radius) : radius(_radius)
{
}

DeltaVelocity::DeltaVelocity(float _refVelocity,float _actualVelo, float _timeDelta) :  referenceVelocity(_refVelocity), actualVelocity(_actualVelo), timeDelta(_timeDelta)
{
}
