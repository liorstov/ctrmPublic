#include "ImageP.h"

ImageP::ImageP(int _index, float _x, float _y, float _z, int jbeg,int jend, int vrange,int dv,int window)
{
	x = _x;
	y = _y;
	z = _z;   
	IPindex = _index;

	for (int i = jbeg+window; i <= jend-window; i++)
	{

		samples.push_back(sampleData(i, vrange,dv));
	}
}


sampleData::sampleData(int number, int vrange, int dv): sampN(number),semblance(0.0f)
{
	for (int i = -vrange; i <= vrange; i+=dv)
	{
		velo.push_back(velocity(i));
	}
}


velocity::velocity(int _value): value(_value), semb(0.0f)
{
}
