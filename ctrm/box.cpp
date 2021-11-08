

#include "box.h"


const char* buildString = "This build XXXX was compiled at  " __DATE__ ", " __TIME__ ".";


#include <sstream> 
box::box(int traces,int lines,int tracesMin, int linesMin,int startRad, int endRad, int _dx, int _dy,int _jbeg, int _jend,int _vrange, int _dv,int _minDist,int _windowSize, int _dr):xmax(traces),ymax(lines),xmin(tracesMin),ymin(linesMin), windowSize(_windowSize), startRadius(startRad), endRadius(endRad), dyLines(_dy),dxTrace(_dx),vRange(_vrange), dv(_dv),jbeg(_jbeg), jend(_jend),minDist(_minDist),dr(_dr)
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
		<< endl;
	/*if (geom.rows() != energy.cols()) {
		throw(std::range_error("number of geo and traces mismatch"));
	}
	std::cout << (geom.rows(), energy.cols());*/
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
		++serialNumber;
	}
	this->signalPosition = signalPos;
	if (jbeg == -1 || jend == -1) {
		jbeg = signalPosition - 200;
		jend = signalPosition + 200;

		if (jbeg < 0) jbeg = 0;
		if (jend > (nsamp-100)) jend = nsamp-100;
	}
	std::cout << serialNumber <<" signalPOsition: "<< this->highestEnergy<<" "<< signalPosition << " start to end " << jbeg << "-" << jend << endl;
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
		radius += dr;
		std::cout << num.second << " ";
	}
}


void box::createImageSpace()
{
	std::cout << "lines"<< ymax << "traces" << xmax << "radius" << startRadius << "to" << endRadius<< endl;
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
	ProgressBar pBar(int(vImagePoints.size()), 70);
	#pragma omp parallel default(none) shared(pBar,counter)
	{
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

							if (deltaSample >= nsamp || deltaSample <= 0 || surface > minDist) {
								continue;
							}
							else {

								geoEnergy = Geo.U[deltaSample];
							}

							SembSize += fabsf(geoEnergy);
							S += geoEnergy;
							SS += (geoEnergy * geoEnergy);

							if (geoEnergy >= 0.0f)
							{
								plusCounter++;
							}
							else
							{
								minusCounter++;
							}
							totalGeophones++;
						}
						
						
						totalWindowIterations++;

					}
					//calculate the weight of the gate
					SemblaneWeight = SembSize / totalGeophones / this->highestEnergy;
					
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

			if (deltaTimemax >= nsamp || deltaTimemin <= 0 || surface > minDist) {
				continue;
			}
			Ip.timeDeltas.push_back(tuple<int,float,float>(Geo.index,deltaTimemin/dt,deltaTimemax/dt));
		}
	}
	std::cout << "finished calculating time delta only lol" << endl;

}
///// <summary>
///// calculate semblence for each IP for each sample
///// </summary>
//void box::corrolationOnGeo(bool RP)
//{
//	std::cout << "calculating corrolation"<< vImagePoints.size()<< "image points";
//
//	//INIT PROGRESS BAR
//	progresscpp::ProgressBar progressBar(int(vImagePoints.size()), 70);
//	float  timeDiff{ 0 }, S{ 0 }, SS{ 0 }, d{ 0 }, f1{ 0 }, f2{ 0 }, fCorrolation, geoEnergy{ 0 }, totalIterations{ 0 };
//	int  nTotalSamples{ 0 }, totalSize{ int(vImagePoints.size()) }, i{ 0 }, nNumOfGeoParticipate{ 0 }, progress{ 0 };
//	for (auto& Ip : vImagePoints)
//	{		
//			for (int sample = jbeg; sample <= jend; ++sample)
//			{
//				f1 = 0;
//				f2 = 0;
//
//				for (int velocity = 0; velocity < vRange / dv*2 + 1; ++velocity)
//				{
//					SS = 0;
//					S = 0;
//					nTotalSamples = 0;
//					for (auto& IpGeo : Ip.IpGeoCalc)
//					{
//						//averageWindow
//						for (int window = sample - windowSize; window <= sample + windowSize; window++)
//						{
//							// the time for the current speed
//							timeDiff = window + IpGeo.VelocityRangeTimeDelta[velocity].timeDelta;
//
//							//std::cout << (radius, velocity,timeDiff);
//							// check if inside the sample frame
//							if (timeDiff >= nsamp || timeDiff <= 0 || IpGeo.horizontal > minDist) {
//
//								continue;
//							}
//
//
//							// get the energy from the input according to geophone number and time difference
//							//geoEnergy = getGeophoneEnergy(IpGeo.GeoIndex, );
//							geoEnergy = vGeophones[IpGeo.GeoIndex-1].U[lroundf(timeDiff)];
//
//							if (RP) {
//								geoEnergy /= IpGeo.RadPatternPwave;
//							}
//							// samblence
//							S += geoEnergy;
//							SS += powf(geoEnergy, 2);
//							totalIterations++;
//
//						}
//						nTotalSamples++;
//
//					}
//
//					// calculate samblence for each velocity
//					if (SS != 0)
//					{
//						d = powf(S, 2) / (SS * nTotalSamples);
//					}
//					else d = 0;
//
//					f1 += d * expf(60 * d);
//					f2 += expf(60 * d);
//					
//				}
//
//				//for each sample calculate the corrolation
//				fCorrolation = f1 / f2;
//
//				//for each sample add the radius and correlation
//				//Ip.SampleSemblance.push_back(std::tuple<int,int, float>(sample, nTotalSamples,fCorrolation));
//			}
//			//std::cout << (Ip.IPindex) << endl;
//			// print status bar
//			++progressBar;
//			if (i % int(vImagePoints.size()*0.1) == 0){
//				progressBar.display();
//			}
//			i++;
//	}
//	progressBar.done();
//
//}
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
//
//Eigen::MatrixXd box::getSample()
//{
//
//	Eigen::MatrixXd ret(vImagePoints.size() * vImagePoints[0].SampleSemblance.size(),7);
//	std::cout << ("writing", vImagePoints.size() * vImagePoints[0].SampleSemblance.size(), "records");
//
//	auto i{ 0 }, stamp{ 0 };
//	std::pair<int, float> semblance;
//	for (auto const& IP: vImagePoints)
//	{
//		stamp = 0;
//		for (auto const& samb : IP.SampleSemblance)
//		{		
//			ret(i, 0) = IP.IPindex;
//			ret(i, 1) = IP.x;
//			ret(i, 2) = IP.y;
//			ret(i, 3) = IP.z;
//			ret(i, 4) = std::get<0>(samb);
//			ret(i, 5) = std::get<1>(samb);
//			ret(i, 6) = std::get<2>(samb);
//			++i;
//		}		
//	}
//
//	std::cout << (i);
//	return ret;
//}
//
//Eigen::MatrixXd box::getIP()
//{
//	Eigen::MatrixXd ret(vImagePoints.size() * vImagePoints[0].IpGeoCalc.size() * vImagePoints[0].IpGeoCalc[0].VelocityRangeTimeDelta.size()*(jend-jbeg), 12);
//	//Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(12);
//	std::cout << ("writing IP");
//	int i{ 0 }, j{ 0 };
//	for (auto& Ip : vImagePoints) {
//		for (auto& GeoIp : Ip.IpGeoCalc) {
//			for (auto& velo : GeoIp.VelocityRangeTimeDelta) {
//				for (int sample = jbeg; sample < jend; ++sample){
//					ret(i, j++) = Ip.x;
//					ret(i, j++) = Ip.y;
//					ret(i, j++) = Ip.z;
//					ret(i, j++) = velo.referenceVelocity;
//					ret(i, j++) = GeoIp.GeoIndex;
//					ret(i, j++) = GeoIp.horizontal;
//					ret(i, j++) = GeoIp.vertical;
//					ret(i, j++) = GeoIp.distance;
//					ret(i, j++) = velo.actualVelocity;
//					ret(i, j++) = velo.timeDelta;
//					ret(i, j++) = sample;
//					ret(i, j++) = getGeophoneEnergy(GeoIp.GeoIndex, lroundf(sample+velo.timeDelta))	;
//					i++;
//					j = 0;
//				}
//			}
//		}
//	}
//
//	return ret;
//}


//
//void box::writeIP()
//{
//	std::cout << "writing space to file..." << endl;
//	std::ofstream myfile;
//	myfile.open("ImagePoints.csv");
//	myfile << "Ip.x, Ip.y, Ip.z, GeoIndex, horizontal, vertical, distance, velocity, timedelta" << endl;
//	for (auto& Ip: vImagePoints) {
//		for (auto& GeoIp : Ip.IpGeoCalc) {
//			for (auto& velo : GeoIp.VelocityRangeTimeDelta) {
//				myfile << Ip.x << "," << Ip.y << "," << Ip.z << "," << GeoIp.writeInfo(int(Ip.x)) << "," <<
//					velo.actualVelocity << "," << velo.timeDelta << endl;
//
//			}
//
//		}
//	}
//	myfile.close();
//
//}
//

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
