#include "ImageP.h"

ImageP::ImageP(int _index, float _x, float _y, float _z)
{
	x = _x;
	y = _y;
	z = _z;
	IPindex = _index;

	IpGeoCalc.reserve(sizeof(IpToGeoGeometry) * 100);
	SampleSemblance.reserve(sizeof(float) * 500);
}
