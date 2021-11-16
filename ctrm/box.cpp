

#include "box.h"

const char* buildString = "This build XXXX was compiled at  " __DATE__ ", " __TIME__ ".";


#include <sstream> 
box::box(int traces,int lines,int tracesMin, int linesMin,int startRad, int endRad, int _dx, int _dy,int _jbeg, int _jend,int _vrange, int _dv,int _minDist,int _windowSize, int _dr,float _numOfThreadsPercent,float _sampleRate):xmax(traces),ymax(lines),xmin(tracesMin),ymin(linesMin), windowSize(_windowSize), startRadius(startRad), endRadius(endRad), dyLines(_dy),dxTrace(_dx),vRange(_vrange), dv(_dv),jbeg(_jbeg), jend(_jend),minDist(_minDist),dr(_dr), numOfThreadsPercent(_numOfThreadsPercent), dt(_sampleRate)
{
	vImagePoints.reserve(sizeof(ImageP) * 22000); 
	std::cout << buildString << "xmax " << xmax<<endl

		<< "ymax " << ymax << endl
		<< "ymin " << ymin << endl
		<< "xmin " << xmin << endl
		<< "startRad " << startRad << endl
		<< "endRad " << endRad << endl
		<<" dr "<< dr << endl
		<< " _dx " << _dx << endl
		<< "_dy " << _dy << endl
		<< "minDist " << minDist << endl
		<< "jbeg " << jbeg << endl
		<< "jend " << jend << endl
		<< "_vrange " << _vrange << endl
		<< "_dv " << _dv << endl
		<< "windowsSize " << windowSize << endl
		<< "threadPercentage " << numOfThreadsPercent << endl
		<< "dt " << dt << endl
		<< endl;	
}
void box::readCoordnumpy(string file){
	vector<unsigned long> shape;
	bool fortran_order;
	vector<double> data;
	npy::LoadArrayFromNumpy(file, shape, fortran_order, data);
	std::cout << "coord file open "<< file << endl;
	int signalPos{ 0 },currentGeo{0}, counter{0}, i{0};
	while(counter < shape.at(0)){
		Geophone newGeo(counter++, data[i], data[i+1], data[i+2]);
		//std::cout <<shape.at(0)<<counter<< data[i]<<data[i+1]<< data[i+2]<< endl;
		vGeophones.push_back(newGeo);
		i+=3;
		if (data[i+2] < minimumDepth)
			minimumDepth = float(data[i+2]);

	}
	std::cout << "number of geophones read: "<<  vGeophones.size()<< endl;

}
//
void box::readCoord(string file)
{
	ifstream coordinats;
	coordinats.open(file);
	if (coordinats.is_open()) {
		std::cout << ("coordinate file open ") << file << endl;
	}
	else {
		std::cout << "can't open coordinate file "<< file<< endl;
	}
	string line;
	float index,x, y, z;
	while (std::getline(coordinats, line)) {


		stringstream ss(line);
		string num;
		ss >> index >> x >> y >> z;
		//std::cout << (index,x,y,z);

		Geophone newGeo(index, x, y, z);
		vGeophones.push_back(newGeo);

		if (z < minimumDepth)
			minimumDepth = float(z);
	}

	std::cout << "number of geophones read: "<<  vGeophones.size()<< endl;
}


void box::readEnergy(string file)
{
	ifstream Ener;
	Ener.open(file,ios::binary);
	if (Ener.is_open()) {
		std::cout << ("Ener file open ") << file << endl;
	}
	else {
		std::cout << "can't open Ener file" << file << endl;
		return;
	}
	int serialNumber{ 0 }, currentGeo{ 0 }, currentSemp{ 0 }, nofGeo{ static_cast<int>(vGeophones.size()) }, signalPos{ 0 };
	
	float value;

	while (Ener.read(reinterpret_cast<char*>(&value), sizeof(value)) && currentGeo < vGeophones.size())
	{
		vGeophones[currentGeo].U.push_back((value));
		value = fabsf(value);
		currentGeo = serialNumber / nsamp;
		currentSemp = serialNumber % nsamp;
		if (value > this->highestEnergy) {
			this->highestEnergy = value;
			signalPos = currentSemp;
		}
		std::cout << currentGeo << " "<< currentSemp << " " <<  value << endl;
		++serialNumber;
	}
	this->signalPosition = signalPos;
	if (jbeg == -1 || jend == -1) {
		jbeg = signalPosition - 200;
		jend = signalPosition + 200;

		if (jbeg < 0) jbeg = 0;
		if (jend > (nsamp-100)) jend = nsamp-100;
	}
	std::cout << "Number of recievers: "<<currentGeo<< "Number of samples: "<<currentSemp <<" high signal value: "<< this->highestEnergy<<" high signal position "<< signalPosition << " start to end " << jbeg << "-" << jend << endl;
}
void box::readEnergyFromNpy(string file)
{
	 vector<unsigned long> shape;
  bool fortran_order;
  vector<double> data;
   npy::LoadArrayFromNumpy(file, shape, fortran_order, data);
   
   
  
   int signalPos{ 0 },currentGeo{0}, counter{0}, samples{int(shape.at(1))};
   for (float value: data)
   {
	   currentGeo = counter / samples;
	   vGeophones[currentGeo].U.push_back((value));
	   if (value > this->highestEnergy) {
			this->highestEnergy = value;
			signalPos = 0;
		}
	  counter++;
	  	//  std::cout << "Number of recievers: "<<currentGeo << "   "<<counter % samples<< " " << value<< endl;
		//   << "Number of samples: "<<currentSemp <<" high signal value: "<< this->highestEnergy<<" high signal position "<< signalPosition << " start to end " << jbeg << "-" << jend << endl;

	}

	std::cout <<endl<< "Number of recievers: "<<currentGeo << endl;
	for (float i: shape)
    std::cout << i << ' ';

}
void box::readVelo(string file)
{
	ifstream velo;
	velo.open(file);
	if (velo.is_open()) {
		std::cout << ("velo file open ") << file << endl;
	}
	else {
		std::cout << "can't open velo file  set to constant" << file << endl;
		for (size_t i = 1; i < 41; i++)
		{
			std::pair<float, float> num{ i,2000 };
			vRadiusVelo.push_back(num);
			std::cout << to_string(vRadiusVelo[i].first) << endl;
		}
		return;
	}
	std::pair<float, float> num;
	float radius = 0;
	while (velo >> num.second)
	{
		num.first = radius;
		vRadiusVelo.push_back(num);
		radius += 1;
		std::cout << num.second << " ";
	}
	if (vRadiusVelo.size() < (endRadius - startRadius)) {
		char buffer [402];
		sprintf(buffer, "number of velocities in velo profile (%d) is lower then cube z dimension (%d)", int(vRadiusVelo.size()), endRadius-startRadius);
		throw(buffer);
	}
}


void box::createImageSpace()
{
	std::cout << "lines: "<< ymax << "traces: " << xmax << "radius: " << startRadius << "to: " << endRadius<< endl;
	int index{ 0 };
	for (float y_index = ymin; y_index <= ymax ; y_index+= dyLines)
	{
		for (float x_index = xmin; x_index <= xmax; x_index+= dxTrace)
		{
			for (float z_index = startRadius; z_index <= endRadius; z_index += dr)
			{
				ImageP newIP(index++, (x_index), ((y_index)), ((z_index)),jbeg,jend,vRange,dv,windowSize);
				vImagePoints.push_back(newIP);
			}
		}
	}
	std::cout << "created: " << vImagePoints.size() << " points" << endl;
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
	float retValue {0};
	auto samples = vGeophones.at(size_t(index)).U;
	if ((timeDiff >= 0) && (timeDiff<samples.size())){
		retValue = samples.at(size_t(timeDiff));
	}
	return retValue;
}

float box::getGeoZ(int x, int y) {
	for (auto const& value : vGeophones) {
		if (value.x == x && value.y == y) {
			return float(value.z);
		}
	}
	return 0;
}

/// <summary>
/// calculates the distances between geophones and Ip 
/// for each pair calculates the potential time differences for each possible velocity
/// </summary>
void box::CalcSurfaceDist() 
{
	int counter{ 0 };
	std::cout << "calculating Ip to Geophone distances" << endl;
	std::cout << "thread percentage: "<<numOfThreadsPercent*100 << endl;
	ProgressBar pBar(int(vImagePoints.size()), 70);

	#pragma omp parallel default(none) shared(pBar,counter) num_threads(int(omp_get_max_threads()*numOfThreadsPercent))
	{
		printf("num of threads %d \n", omp_get_num_threads( ));
		float  distance{ 0 }, surface{ 0 }, VVa{ 0 }, VV, CurrectVelocity, IpDepth{ 0 }, ydist{ 0 }, xdist{ 0 }, geoEnergy{ 0 };
		float  minusCounter{ 1 }, plusCounter{ 1 },SembSize{ 0 },totalGeophones{ 0 }, SemblaneWeight{ 1 }, deltaTime{ 0 }, S{ 0 }, SS{ 0 }, f1{ 0 }, f2{ 0 }, fCorrolation{ 0 }, windowAvr{ 0 }, currentSemblance{ 0 };
		int   deltaSample{ 0 },totalSize{ int(vImagePoints.size()) }, totalWindowIterations{ 0 };
		auto GeoVector = vGeophones;

		#pragma omp for  
		for (int i = 0; i < vImagePoints.size(); i++)
		{

			auto& Ip = vImagePoints[i];
			fCorrolation = 0;
			for (auto& sample : Ip.samples)
			{
				f1 = f2 = 0;
				for (auto& velocity : sample.velo)
				{

					totalWindowIterations = 0;
					windowAvr = 0;
					totalGeophones = S = SS = SembSize = 0.0f;
					SemblaneWeight = 1;
					minusCounter = plusCounter = 0.0f;
					for (int window = sample.sampN - windowSize; window <=  sample.sampN+ windowSize; window++)
					{
						//loop over geophones
						for (auto& Geo : GeoVector)
						{
							//geophone depth
							Geo.z;
							IpDepth = (Ip.z - Geo.z);
							ydist = Ip.y - Geo.y;
							xdist = Ip.x - Geo.x;
							surface = sqrtf(powf(Ip.x - Geo.x, 2) + powf(Ip.y - Geo.y, 2));
							distance = sqrtf(powf(surface, 2) + powf(IpDepth, 2));

							//distance between Ip and Geo for current radius
							CurrectVelocity = vRadiusVelo.at(int(Ip.z)).second;
							//average velocity from the Ip to the Geo
							VVa = calcAvarageVelo(Ip.z, Geo.z);
							// a loop for checking a range of values around the speed of the radius 

							// the avarege value plus an offset
							VV = VVa + velocity.value;
							deltaTime = distance / VV - float(Ip.z) / CurrectVelocity;
							deltaTime = deltaTime / dt;

							deltaSample = window + int(roundf(deltaTime));

							if (deltaSample <= 0 || surface > minDist) {
								continue;
							}
							else {

								geoEnergy = Geo.U[deltaSample];
							}

							SembSize += fabsf(geoEnergy);
							S += geoEnergy;
							SS += (geoEnergy * geoEnergy);

							totalGeophones++;
						}
						totalWindowIterations++;
					}
					//calculate the weight of the gate
					//SemblaneWeight = SembSize / totalGeophones / this->highestEnergy;
					
					//number of geophones is calculated again not including the window
					totalGeophones = totalGeophones/ totalWindowIterations;
					
					//calculate semblance for the depth
					if (SS != 0)						
							currentSemblance = (S*S) / (SS * totalGeophones);
						else
							currentSemblance = 0.0f;

					velocity.semb = currentSemblance* SemblaneWeight;
					/*if (minusCounter>0)
						minusCounter = minusCounter / totalGeophones * logf(minusCounter / totalGeophones);
					if (plusCounter>0)
						plusCounter = plusCounter / totalGeophones * logf(plusCounter / totalGeophones);
					velocity.semb =  (1-(plusCounter+minusCounter))*SemblaneWeight;*/
					f1 += velocity.semb * expf(60 * velocity.semb);
					f2 += expf(60 * velocity.semb);
					/*if (Ip.x == 30 && Ip.y == -25 && Ip.z == 1) {
						std::cout << velocity.value << " " << f1 << " " << f2 << " " << velocity.semb << endl;
					}*/
				}
					
				sample.semblance = f1/f2;
			}
			counter++;
			++pBar;
			if (counter % int(totalSize * 0.01) == 0) {
				//std::cout<< counter << " out of "<< totalSize<<  endl;
				#pragma  omp critical
				pBar.display();
			}
		}
		
	}
	pBar.done();
	std::cout << "done" << endl;
}
void box::CalcTimeDeltaOnly()
{
	std::cout << "calculating time delta only lol" << endl;

	float  distance{ 0 }, surface{ 0 }, VVmax, VVmin, CurrectVelocity, IpDepth{ 0 }, ydist{ 0 }, xdist{ 0 }, geoEnergy{ 0 };
	float  VV{ 0 }, S{ 0 }, SS{ 0 }, f1{ 0 }, f2{ 0 }, fCorrolation{ 0 }, windowAvr{ 0 }, currentSemblance{ 0 };
	int  totalGeophones{ 0 }, totalSize{ int(vImagePoints.size()) }, totalWindowIterations{ 0 };
	for (auto& Ip : vImagePoints) {
		for (auto& Geo : vGeophones)
		{
			auto deltaTimemin = 0.0f;
			auto deltaTimemax = 0.0f;
			//geophone depth
			Geo.z;
			IpDepth = (Ip.z - Geo.z);
			ydist = Ip.y - Geo.y;
			xdist = Ip.x - Geo.x;
			surface = sqrtf(powf(Ip.x - Geo.x, 2) + powf(Ip.y - Geo.y, 2));
			distance = sqrtf(powf(surface, 2) + powf(IpDepth, 2));


		

			//distance between Ip and Geo for current radius
			CurrectVelocity = vRadiusVelo.at(int(Ip.z)).second;
			//average velocity from the Ip to the Geo
			VV = calcAvarageVelo(Ip.z, Geo.z);
			// a loop for checking a range of values around the speed of the radius 

			VVmax = VV+ vRange;
			VVmin = VV-vRange;
			deltaTimemax = distance / VVmax - float(Ip.z) / CurrectVelocity;
			deltaTimemin = distance / VVmin - float(Ip.z) / CurrectVelocity;

			if (deltaTimemin <= 0 || surface > minDist) {
				continue;
			}
			Ip.timeDeltas.push_back(tuple<int,float,float>(Geo.index,deltaTimemin/dt,deltaTimemax/dt));
		}
	}
	std::cout << "finished calculating time delta only lol" << endl;

}

float box::calcAvarageVelo(float IpDepth, float GeoDepth)
{
	auto begin = min(int(roundf(IpDepth)), int(roundf(GeoDepth)));
	auto end = max(int(roundf(IpDepth)), int(roundf(GeoDepth)));
	float  CurrentVelocity, avarage = 0, counter{ 0 };
	int radius;
	for (auto const& value : vRadiusVelo)
	{
		radius = value.first;
		CurrentVelocity = value.second;

		if (radius >= begin && radius <= end) {
			avarage += (1 / CurrentVelocity);
			counter++;
		}
		
	}

	return((counter) / avarage);
}

void box::writeSemblenceNpy(std::string file)
{
	std::cout << "creating output file cube" << file << endl;

	auto i{ 0 };
	vector<float> ret;
	ret.reserve(vImagePoints.size() * vImagePoints[0].samples.size() * 5);
	const long unsigned leshape[] = { vImagePoints.size() * vImagePoints[0].samples.size(),5 };

	for (auto& Ip : vImagePoints) {
		for (auto& sample : Ip.samples) {
			ret.push_back(Ip.x);
			ret.push_back(Ip.y);
			ret.push_back(Ip.z);
			ret.push_back(sample.sampN);
			ret.push_back(sample.semblance);
		}
	}
	std::cout << "save to file: " << file << endl;
	npy::SaveArrayAsNumpy(file, false, 2, leshape, ret);

}

void box::writeTimeDeltasNpy(std::string file)
{
	std::cout << "creating output file delta" << file << endl;

	auto i{ 0 };
	vector<float> ret;
	ret.reserve(vImagePoints.size() * vImagePoints[0].timeDeltas.size() * 6);
	const long unsigned leshape[] = { vImagePoints.size() * vImagePoints[0].timeDeltas.size(),6 };

	for (auto& Ip : vImagePoints) {
		for (auto& Delta : Ip.timeDeltas) {
			ret.push_back(std::get<0>(Delta));
			ret.push_back(Ip.x);
			ret.push_back(Ip.y);
			ret.push_back(Ip.z);
			ret.push_back(std::get<1>(Delta));
			ret.push_back(std::get<2>(Delta));
		}
	}
	npy::SaveArrayAsNumpy(file, false, 2, leshape, ret);
}
