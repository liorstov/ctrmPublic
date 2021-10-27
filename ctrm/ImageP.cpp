#include "ImageP.h"

ImageP::ImageP(int _index, float _x, float _y, float _z, int jbeg,int jend, int vrange,int dv,int window)
{
	x = _x;
	y = _y;
	z = _z;   
	IPindex = _index;

	for (int i = jbeg+window; i <= jend-window; i++)
	{

		samples.push_back(sample(i, vrange,dv));
	}
}

char* ImageP::toChar()
{
	char * array = new char[sizeof(float) * (3+ samples.size())];
	char buffer[sizeof(float)*3];
	sprintf(buffer, "%.4f%.4f%.4f", x,y,z);
	sprintf(buffer, "%f", y);
	sprintf(buffer, "%f", z);

	for(auto &sample : this->samples)
	{
		sprintf(array, sample.toChar());
	}
	return array;
}

sample::sample(int number, int vrange, int dv): sampN(number),semblance(0.0f)
{
	for (int i = -vrange; i <= vrange; i+=dv)
	{
		velo.push_back(velocity(i));
	}
}

char* sample::toChar()
{
	char array[sizeof(int) + sizeof(float)];
	sprintf(array, "%d", this->sampN);
	sprintf(array, "%f", this->semblance);
	return array;
}

velocity::velocity(int _value): value(_value), semb(0.0f)
{
}
