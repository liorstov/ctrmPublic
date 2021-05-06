#include <pybind11/pybind11.h>

#include "box.h"
namespace py = pybind11;


#include <sstream> 
box::box(int lines,int traces,int startRad, int endRad):numlin(lines), numcen(traces), startRadius(startRad), endRadius(endRad)
{
	vGeophones.reserve(sizeof(Geophone) * 10000);
	vImagePoints.reserve(sizeof(ImageP) * 10000);

	
}
void box::readCoord(string file)
{
	py::print("reading lol");
	ifstream coordinats;
	coordinats.open("OneHorzShot4CTRM_R_Coord.txt");
	if (coordinats.is_open()) {
		py::print("file open");
	}
	else {
		py::print("no open");
	}
	string line;
	int index,x, y, z;
	for (size_t i = 0; i < numrecl; i++)
	{
		for (size_t j = 0; j < ntrace; j++)
		{
			getline(coordinats, line);

			stringstream ss(line);
			string num;
			ss >> index >> x >> y >> z;
			//py::print(index,x,y,z);

			Geophone newGeo(index, x, y, z);
			vGeophones.push_back(newGeo);

			if (z < minimumDepth)
				minimumDepth = float(z);
		}

	} 
}

void box::readVelo(string file)
{
	ifstream velo;
	velo.open("vrs.txt");
	std::pair<float,float> num;
	float radius =rmin;
	while (velo >> num.second)
	{
		num.first = radius;
		vRadiusVelo.push_back(num);
		radius += dr;
	}
}

void box::readEnergy(string file)
{
	ifstream Ener;
	Ener.open("Zc.nor",ios::binary);
	int serialNumber{ 0 }, currentGeo{ 0 };
	float value;
	while (Ener.read(reinterpret_cast<char*>(&value), sizeof(value)) && currentGeo < vGeophones.size())
	{
		vGeophones[currentGeo].U.push_back((value));
		currentGeo = ++serialNumber / nsamp;
	}
}

void box::createImageSpace()
{
	py::print("lines", numlin, "traces", numcen, "radius", startRadius, "to", endRadius);
	int index{ 0 };
	for (int y_index = 0; y_index < numcen; y_index++)
	{
		for (int x_index = 0; x_index < numlin; x_index++)
		{
			ImageP newIP(index++, float(x_index+1) * dxc, (float(y_index))* dyc + coordcy0, getGeoZ(x_index,y_index));
			vImagePoints.push_back(newIP);
		}
	}

}

float box::getIP(int point, int sample)
{
	return std::get<0>(vImagePoints[point].SampleSemblance[sample]);
}

std::vector<int> box::getCoord()
{
	std::vector<int> ret;
	for (auto const& value : vGeophones) {
		ret.push_back(value.index);
	}
	return ret;
}

float box::getGeophoneEnergy(int index, int timeDiff)
{	
	return vGeophones[size_t(index) - 1].U[size_t(timeDiff)-1];
}

float box::getGeoZ(int x, int y) {
	for (auto const& value : vGeophones) {
		if (value.x == x && value.y == y) {
			return float(value.z);
		}
	}
	return 0;
}

void box::CalcSurfaceDist() 
{
	py::print("calculating Ip to Geophone distances");

	float dist, betta, verDepth, VVa{ 0 }, VV, radius, CurrectVelocity, IpDepth{ 0 };
	int deltaTime{ 0 };
	for (auto& Ip : vImagePoints)
	{
		auto start = std::chrono::high_resolution_clock::now();

		for (auto const& Geo : vGeophones)
		{
			dist = sqrtf(powf(Ip.x - Geo.x, 2) + powf(Ip.y - Geo.y, 2) + powf(Ip.z - Geo.z, 2));
			verDepth = Ip.z - minimumDepth + Geo.z;
			betta = asinf(verDepth / dist);
			auto newIpGeo = IpToGeoGeometry(Geo.index, dist, betta, verDepth);

			for (auto const& value : vRadiusVelo) {	
				radius = value.first;
				CurrectVelocity = value.second;
				if (radius < startRadius || radius > endRadius) continue;			
				

				auto newRadiusVelo = RadiusTimeDelta(radius);
				VVa = calcAvarageVelo(verDepth, radius);
				IpDepth = verDepth + radius;
				//distance between Ip and Geo for current radius
				newRadiusVelo.distanceToIP = sqrtf(powf(dist, 2) + powf(IpDepth, 2) - 2 * sinf(betta) * dist * (IpDepth));
				// a loop for checking a range of values around the speed of the radius 
				for (size_t i = 0; i < vRange/dv*2+1; i++)
				{
					// the avarege value plus an offset
					VV = VVa + i * dv - vRange;
					deltaTime = lroundf((newRadiusVelo.distanceToIP / VV - (IpDepth) / CurrectVelocity) / dt);
					newRadiusVelo.VelocityRangeTimeDelta.push_back(DeltaVelocity(VV,deltaTime));
				}
				newIpGeo.RadiusTimeValues.push_back(newRadiusVelo);
			}			
			Ip.IpGeoCalc.push_back(newIpGeo);	

		}
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		//std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	}

	
}
void box::corrolationOnGeo()
{
	cout << "calculating corrolation" << endl;
	float  S{ 0 }, SS{ 0 }, d{ 0 }, f1{ 0 }, f2{ 0 }, fCorrolation, geoEnergy{ 0 }, progress{ 0 }, totalIterations{ 0 };
	int timeDiff, nTotalTraces{ 0 }, totalSize{ int(vImagePoints.size()) }, i{ 0 };

	for (auto& Ip : vImagePoints)
	{
		for (int radius = 0; radius <= endRadius-startRadius; radius++)
		{
			for (int sample = jbeg; sample < jend; sample += lroundf(resamp))
			{
				f1 = 0;
				f2 = 0;

				for (int velocity = 0; velocity < vRange / dv*2 + 1; velocity++)
				{
					SS = 0;
					S = 0;
					nTotalTraces = 0;
					for (auto& IpGeo : Ip.IpGeoCalc)
					{
						// the time for the current speed
						timeDiff = sample + IpGeo.RadiusTimeValues[radius].VelocityRangeTimeDelta[velocity].timeDelta;

						// check if inside the sample frame
						if (timeDiff >= nsamp || timeDiff <= 0 || IpGeo.GeoIndex > numrecl * ntrace || IpGeo.surfaceDist > offmax) 
							continue;
						
						// get the energy from the input according to geophone number and time difference
						geoEnergy = getGeophoneEnergy(IpGeo.GeoIndex, int(timeDiff));
						// samblence
						S += geoEnergy;
						SS += powf(geoEnergy, 2);
						nTotalTraces++;
						totalIterations++;
					}

					// calculate samblence for each velocity
					if (SS != 0)
					{
						d = powf(S, 2) / (SS * nTotalTraces);
					}
					else d = 0;

					f1 += d * expf(60 * d);
					f2 += expf(60 * d);

				}

				//for each sample calculate the corrolation
				fCorrolation = f1 / f2;
				Ip.SampleSemblance.push_back(std::tuple<int,int, float>(sample,radius+startRadius,fCorrolation));
			}

		}
		progress = (float)i++ / (float)totalSize;
		if (fmodf(progress, 0.1f) == 0) {
			char str[80];
			std::sprintf(str, "%d%% ", int(floorf(progress * 100)));
			py::print(str);
		}
		
	}

	py::print("Calculated", totalIterations,"iterations");

}
float box::calcAvarageVelo(float begin, float end)
{
	float radius, CurrectVelocity, avarage = 0, counter{ 0 };
	for (auto const& value : vRadiusVelo)
	{
		radius = value.first;
		CurrectVelocity = value.second;

		if (radius > begin && radius <= end) {
			avarage += (1 / CurrectVelocity);
			counter++;
		}
	}

	return((counter) / avarage);
}

Eigen::MatrixXd box::getSample()
{

	Eigen::MatrixXd ret(numcen*numlin*(endRadius-startRadius+1)*nsamp,5);
	
	auto i{ 0 }, stamp{ 0 };
	std::pair<int, float> semblance;
	for (auto const& IP: vImagePoints)
	{
		stamp = 0;
		for (auto const& samb : IP.SampleSemblance)
		{		
			ret(i, 0) = IP.x;
			ret(i, 1) = IP.y;
			ret(i, 2) = std::get<0>(samb);
			ret(i, 3) = std::get<1>(samb);
			ret(i, 4) = std::get<2>(samb);
			++i;
		}		
	}

	return ret;
}

void box::writeIP()
{
	cout << "writing space to file..." << endl;
	std::ofstream myfile;
	myfile.open("ImagePoints.csv");
	myfile << "Ip.x, Ip.y, Ip.z, GeoIndex, surfaceDist, bettaAngle, verDepth, radius,IpGeoDist, velocity, timedelta" << endl;
	for (auto& Ip: vImagePoints) {
		for (auto& GeoIp : Ip.IpGeoCalc) {
			for (auto& radius : GeoIp.RadiusTimeValues) {
				for (auto& velo : radius.VelocityRangeTimeDelta) {

					myfile << Ip.x << "," << Ip.y << "," << Ip.z << "," << GeoIp.writeInfo(int(Ip.x)) << ","<<radius.radius << "," << radius.distanceToIP << "," <<
						velo.actualVelocity << "," << velo.timeDelta << endl;
				}
			}

		}
	}
	myfile.close();

}

void box::writeSemblence()
{
	int i{ 0 };
	std::ofstream myfile;
	myfile.open("samb.csv");
	myfile << "Ip.x, Ip.y, Ip.z, sample number,radius, semblance"<< endl;
	for (auto& Ip : vImagePoints) {
		i = 0;
		for (auto& sample : Ip.SampleSemblance) {
			myfile << Ip.x << "," << Ip.y << "," << Ip.z << "," << std::get<0>(sample)<< ","<< std::get<1>(sample) << "," << std::get<2>(sample) << endl;;
		}
	}
	myfile.close();
}

