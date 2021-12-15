#include "box.h"
#include <sstream> 

const char* buildString = "This build XXXX was compiled at  " __DATE__ ", " __TIME__ ".";

// constructor of the box takes the following parameters:
// xmin,ymin,xmax,ymax dimension of cube
// endrad startrad: depth of cube
//dx,dy,dz cube resolution
// minDist: all recievers within this range should be calculated for each imagePoint
//window size: the size of the semblance gate
//jbeg jend the sesmogram time slice for computation
//numOfThreadsPercent: the fraction of threads to be used
//_sampleRate sesmogram sample rate in miliseconds (seconds/samples)

box::box(int _xmax,int _ymax,int _xmin, int _ymin,int startRad, int endRad, int _dx, int _dy,int _jbeg, int _jend,int _vrange, int _dv,int _minDist,int _windowSize, int _dr,float _numOfThreadsPercent,float _sampleRate):xmax(_xmax),ymax(_ymax),xmin(_xmin),ymin(_ymin), windowSize(_windowSize), startRadius(startRad), endRadius(endRad), dyLines(_dy),dxTrace(_dx),vRange(_vrange), dv(_dv),jbeg(_jbeg), jend(_jend),minDist(_minDist),dr(_dr), numOfThreadsPercent(_numOfThreadsPercent), dt(_sampleRate)
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

//read coordinates of geophones from a numpy array of index,x,y,z for each reciever
void box::readCoordnumpy(string file){

	// read numpy
	vector<unsigned long> shape;
	bool fortran_order;
	vector<double> data;
	std::cout << "opening coord file "<< file << endl;
	npy::LoadArrayFromNumpy(file, shape, fortran_order, data);	
	std::cout << "coord file open "<< file << endl;
	
	int signalPos{ 0 },currentGeo{0}, counter{0}, i{0};

	// move to geophone vector
	while(counter < shape.at(0)){
		Geophone newGeo(data[i], data[i+1], data[i+2], data[i+3]);
		counter++;
		//std::cout <<shape.at(0)<<counter<< data[i]<<data[i+1]<< data[i+2]<< endl;
		vGeophones.push_back(newGeo);
		i+=4;
		if (data[i+3] < minimumDepth)
			minimumDepth = float(data[i+3]);

	}
	std::cout << "number of geophones read: "<<  vGeophones.size()<< endl;

}

// read sesmogram from a numpy array in segy format [tracesXsamples]
void box::readEnergyFromNpy(string file)
{
	//read file
	vector<unsigned long> shape;
	bool fortran_order;
	vector<float> data;
	cout << "read sesmogram from: " << file << endl;
	npy::LoadArrayFromNumpy(file, shape, fortran_order, data);   
	cout << "shape of sesmogram is: " << shape[0] << "  " << shape[1] << endl;

   int signalPos{ 0 },currentGeo{0}, counter{0}, samples{int(shape.at(1))},currentSemp{ 0 };
   nsamp = samples;
   for (auto value: data)
   {
		currentGeo = counter / samples;

		// read data for each  Geophone
		vGeophones[currentGeo].U.push_back((value));
		if (value > this->highestEnergy) {
				// save the highest energy read
				this->highestEnergy = value;
				signalPos = 0;
			}
		counter++;
	}

	// read all if -1 
	if (jbeg == -1 || jend == -1) {
		jbeg = 0;
		jend = samples;
	}

	std::cout <<endl<< "Number of recievers: "<<vGeophones.size() << endl<< "highest energy: " << this-> highestEnergy << endl;;

	// check if number of traces equal number of geophones
	if (vGeophones.size() != currentGeo+1){
		throw(runtime_error("warning : size of sesmogram different from number of recievers"));
	}
	for (float i: shape)
    std::cout << i << ' ';
}

// read velocity from a numpy array of one dimension: velocity for each meter
void box::readVelo(string file)
{
	vector<unsigned long> shape;
	bool fortran_order;
	vector<double> data;
	npy::LoadArrayFromNumpy(file, shape, fortran_order, data); 
	cout << "read velocity profile from "<< file << endl;
	float radius = 0;
	for (auto velo : data){
		std::pair<float, float> num;
		num.first = radius++;
		num.second = velo;
		vRadiusVelo.push_back(num);
		std::cout << velo<< endl;
	}
	std::cout << "finished"<< endl;

	// check if velocity profile fit each computed depth
	if (vRadiusVelo.size() < (endRadius - startRadius)) {
		char buffer [402];
		sprintf(buffer, "number of velocities in velo profile (%d) is lower then cube z dimension (%d)", int(vRadiusVelo.size()), endRadius-startRadius);
		throw(buffer);
	}
}

// create vector of image space
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
		float  minusCounter{ 1 }, plusCounter{ 1 },SembSize{ 0 },totalGeophones{ 0 },totalGeophonesWeight{ 0 }, SemblaneWeight{ 1 }, deltaTime{ 0 }, S{ 0 }, SS{ 0 }, f1{ 0 }, f2{ 0 }, fCorrolation{ 0 }, windowAvr{ 0 }, currentSemblance{ 0 };
		int   deltaSample{ 0 },totalSize{ int(vImagePoints.size()) }, totalWindowIterations{ 0 };
		auto GeoVector = vGeophones;

		// prallel computation of the image space
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
					totalGeophonesWeight = 0;
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

							if (deltaSample >= nsamp || deltaSample <= 0 ) {
								geoEnergy = 0;
							}
							else {
								geoEnergy = Geo.U[deltaSample];
							}

							if (surface<50){
								SembSize += (fabsf(geoEnergy));
								totalGeophonesWeight++;
							}
							S += geoEnergy;
							SS += (geoEnergy * geoEnergy);

							totalGeophones++;
						}
						totalWindowIterations++;
					}

					//number of geophones is calculated again not including the window
					totalGeophones = totalGeophones/ totalWindowIterations;

					//calculate the weight of the gate
					SemblaneWeight = SembSize / totalGeophonesWeight;
					SemblaneWeight = 1;
								
					
					//calculate semblance for the depth
					if (SS != 0)						
							currentSemblance = (S*S) / (SS * totalGeophones);
						else
							currentSemblance = 0.0f;

					velocity.semb = currentSemblance* SemblaneWeight;
					
					// multipass for velocity 
					f1 += velocity.semb * expf(60 * velocity.semb);
					f2 += expf(60 * velocity.semb);
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

// calculatimng only the time deltas for each Ip
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

//Helper function for calculating average velocity for a specific depth
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


//write outpust to numoy
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

// write output of time deltas to seperate file
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
