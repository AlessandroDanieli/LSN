#include <cmath>
#include "stat.h"
using namespace std;

Stat :: Stat(): rnd(nullptr), Nsample(0), Nblocks(1) {}

Stat :: Stat(unsigned int _Nsample, unsigned int _Nblocks) {
	if (_Nsample < _Nblocks) {
		cerr << "Error [Stat :: Stat(unsigned int, unsigned int)]: Nsample < Nblocks" << endl;
		exit(-1);
	}
	else if (_Nsample % _Nblocks != 0) {
		cerr << "Error [Stat :: Stat(unsigned int, unsigned int)]: Nsample/Nblocks is not an integer" << endl;
		exit(-1);
	}
	else if (_Nblocks == 0) {
		cerr << "Error [Stat :: Stat(unsigned int, unsigned int)]: Nblocks = 0" << endl;
		exit(-1);
	}
	else {
		Nsample = _Nsample;
		Nblocks = _Nblocks;
	}
}

Stat :: Stat(Random* r, unsigned int _Nsample, unsigned int _Nblocks) {
	rnd = r;
	if (_Nsample < _Nblocks) {
		cerr << "Error [Stat :: Stat(Random*, unsigned int, unsigned int)]: Nsample < Nblocks" << endl;
		exit(-1);
	}
	else if (_Nsample % _Nblocks != 0) {
		cerr << "Error [Stat :: Stat(Random*, unsigned int, unsigned int)]: Nsample/Nblocks is not an integer" << endl;
		exit(-1);
	}
	else if (_Nblocks == 0) {
		cerr << "Error [Stat :: Stat(Random*, unsigned int, unsigned int)]: Nblocks = 0" << endl;
		exit(-1);
	}
	else {
		Nsample = _Nsample;
		Nblocks = _Nblocks;
	}
}

Stat :: Stat(const vector<double>& _data): rnd(nullptr) {
	for (unsigned int i = 0; i < _data.size(); i++) {
		data.push_back(_data[i]);
	}
	Nsample = _data.size();
	Nblocks = _data.size();
}
	
Stat :: Stat(const vector<double>& _data, unsigned int _Nsample, unsigned int _Nblocks): rnd(nullptr) {
	if (_Nsample < _Nblocks) {
		cerr << "Error [Stat :: Stat(const vector<double>&, unsigned int, unsigned int)]: Nsample < Nblocks" << endl;
		exit(-1);
	}
	else if (_Nsample % _Nblocks != 0) {
		cerr << "Error [Stat :: Stat(const vector<double>&, unsigned int, unsigned int)]: Nsample/Nblocks is not an integer" << endl;
		exit(-1);
	}
	else if (_Nblocks == 0) {
		cerr << "Error [Stat :: Stat(const vector<double>&, unsigned int, unsigned int)]: Nblocks = 0" << endl;
		exit(-1);
	}
	else if (_Nsample > _data.size()) {
		cerr << "Error [Stat :: Stat(const vector<double>&, unsigned int, unsigned int)]: Nsample exceeds data.size() " << endl;
		exit(-1);
	}
	else {
		Nsample = _Nsample;
		Nblocks = _Nblocks;
		for (unsigned int i = 0; i < _data.size(); i++) {
			data.push_back(_data[i]);
		}
	}
}

Stat :: Stat(const Stat& s) {
	rnd = s.GetRandomGen();
	data = s.GetData();
	Nsample = s.GetSampleSize();
	Nblocks = s.GetBlocks();
	exp_mean = s.GetExpMean();
}

Stat :: ~Stat() {}

double Stat :: SetData(const vector<double>& x) { // Sets data and Nsample = data size, block lenght = 1
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample > x.size()) {
		cout << "Warning [Stat :: SetData(vector<double>&)]: since Nsample exceeds the copied vector's size, it is shrunk to fit the new size";
		Nsample = x.size();
		if (Nsample % Nblocks != 0) {
			cout << "Warning [Stat :: SetData(vector<double>&)]: since the (new Nsample)/Nblocks is not an integer, Nblocks is set to 1";
			Nblocks = 1;
		}
	}
	data.clear();
	data = x;
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: CopyData(vector<double> x) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample > x.size()) {
		cout << "Warning [Stat :: CopyData(vector<double>&)]: since Nsample exceeds the copied vector's size, it is shrunk to fit the new size";
		Nsample = x.size();
		if (Nsample % Nblocks != 0) {
			cout << "Warning [Stat :: CopyData(vector<double>&)]: since the (new Nsample)/Nblocks is not an integer, Nblocks is set to 1";
			Nblocks = 1;
		}
	}
	
	data.resize(x.size());
	for (unsigned int i = 0; i < x.size(); i++) data[i] = x[i];
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: Set(const vector<double>& _data, unsigned int _Nsample, unsigned int _Nblocks) {
	auto t1 = chrono::high_resolution_clock::now();
	if (_Nsample < _Nblocks) {
		cerr << "Error [Stat :: Set(const vector<double>&, unsigned int, unsigned int)]: Nsample < Nblocks" << endl;
		exit(-1);
	}
	else if (_Nsample % _Nblocks != 0) {
		cerr << "Error [Stat :: Set(const vector<double>&, unsigned int, unsigned int)]: Nsample/Nblocks is not an integer" << endl;
		exit(-1);
	}
	else if (_Nblocks == 0) {
		cerr << "Error [Stat :: Set(const vector<double>&, unsigned int, unsigned int)]: Nblocks = 0" << endl;
		exit(-1);
	}
	else if (_Nsample > _data.size()) {
		cerr << "Error [Stat :: Set(const vector<double>&, unsigned int, unsigned int)]: Nsample exceeds data.size() " << endl;
		exit(-1);
	}
	else {
		data.clear();
		data = _data;
		Nsample = _Nsample;
		Nblocks = _Nblocks;
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

void Stat :: Append(double x) {
	data.push_back(x); 
	return;	
}

double Stat :: FillUniform() {
	if (rnd == nullptr) {
		cerr << "Error [Stat :: FillUniform()]: random generator rnd is not defined" << endl;
		exit(-1);
	}
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (unsigned int i = 0; i < Nsample; i++) {
		data.push_back(rnd->Uniform());
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
double Stat :: FillUniform(double a, double b) {
	if (rnd == nullptr) {
		cerr << "Error [Stat :: FillUniform(double, double)]: random generator rnd is not defined" << endl;
		exit(-1);
	}
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (unsigned int i = 0; i < Nsample; i++) {
		data.push_back(rnd->Uniform(a, b));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
double Stat :: FillExp(double lambda) {
	if (rnd == nullptr) {
		cerr << "Error [Stat :: FillExp(double)]: random generator rnd is not defined" << endl;
		exit(-1);
	}
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (unsigned int i = 0; i < Nsample; i++) {
		data.push_back(rnd->Exp(lambda));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
double Stat :: FillLorentz(double m, double gamma) {
	if (rnd == nullptr) {
		cerr << "Error [Stat :: FillLorentz(double, double)]: random generator rnd is not defined" << endl;
		exit(-1);
	}
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (unsigned int i = 0; i < Nsample; i++) {
		data.push_back(rnd->Lorentz(m, gamma));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
double Stat :: FillGauss(double m, double sigma) {
	if (rnd == nullptr) {
		cerr << "Error [Stat :: FillGauss(double, double)]: random generator rnd is not defined" << endl;
		exit(-1);
	}
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (unsigned int i = 0; i < Nsample; i++) {
		data.push_back(rnd->Gauss(m, sigma));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


/// ****** ****** ****** ****** ****** ******  Data Analysis - Block Methods ****** ****** ****** ****** ****** ****** ///

double Stat :: Average(double& average) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: Average(double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	else {
		average = 0;	
		for (unsigned int k = 0; k < Nsample; k++) average += data[k];
		average /= (double)Nsample;
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: AverageSigma(double& average, double& sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: AverageSigma(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	else {
		average = 0;
		double average2 = 0;
			
		for (unsigned int k = 0; k < Nsample; k++) {
			average += data[k];
			average2 += data[k]*data[k];
		}
		average /= (double)Nsample;
		average2 /= (double)Nsample;
		sigma = sqrt(average2 - average*average);
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: AverageVariance(double& average, double& variance) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: AverageVariance(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	else {
		average = 0;
		double average2 = 0;
			
		for (unsigned int k = 0; k < Nsample; k++) {
			average += data[k];
			average2 += data[k]*data[k];
		}
		average /= (double)Nsample;
		average2 /= (double)Nsample;
		variance = average2 - average*average;
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Calculates the average of each block 
double Stat :: BlockAverage(vector<double>& average) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockAverage(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		average.clear();
		for (unsigned int n = 0; n < Nblocks; n++) {
			average.push_back(0);
			if (L == 1) { // Blocks of length 1
				average[n] = data[n];
			}
			else {
				for (unsigned int k = 0; k < L; k++) {
					average[n] += data[n*L + k];
				}
				average[n] /= (double)L;
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Quadratic deviation from block average, calculated for each block of data
double Stat :: BlockSigma(vector<double>& sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockSigma(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		sigma.clear();
		
		if (L == 1) {
			cout << "Warning [Stat :: BlockSigma(vector<double>&)]: the standard deviation is automatically set to zero";
			cout << " since each block has length 1." << endl;
			for (unsigned int n = 0; n < Nblocks; n++) sigma.push_back(0);
		}
		else {
			double average, average2;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average = 0;
				average2 = 0;
				sigma.push_back(0);
				
				for (unsigned int k = 0; k < L; k++) {
					average += data[n*L+k];
					average2 += data[n*L+k]*data[n*L+k];
				}
				average /= (double)L;
				average2 /= (double)L;
				sigma[n] = sqrt(average2 - average*average);
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


// Quadratic deviation from block average, calculated for each block of data
double Stat :: BlockVariance(vector<double>& variance) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockVariance(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		variance.clear();
		
		if (L == 1) {
			cout << "Warning [Stat :: BlockVariance(vector<double>&)]: the standard deviation is automatically set to zero";
			cout << " since each block has length 1." << endl;
			for (unsigned int n = 0; n < Nblocks; n++) variance.push_back(0);
		}
		else {
			double average, average2;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average = 0;
				average2 = 0;
				variance.push_back(0);
				
				for (unsigned int k = 0; k < L; k++) {
					average += data[n*L+k];
					average2 += data[n*L+k]*data[n*L+k];
				}
				average /= (double)L;
				average2 /= (double)L;
				variance[n] = average2 - average*average;
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Mean and Standard deviation (using block's mean) of each block
double Stat :: BlockAverageSigma(vector<double>& average, vector<double>& sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockAverageSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		average.clear();
		sigma.clear();
		if (L == 1) {
			cout << "Warning [Stat :: BlockAverageSigma(vector<double>&, vector<double>&)]: the standard deviation is automatically set to zero";
			cout << " since each block has length 1." << endl;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average.push_back(data[n]);
				sigma.push_back(0);
			}
		}
		else {
			double average2;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average.push_back(0);
				average2 = 0;
				sigma.push_back(0);
				
				for (unsigned int k = 0; k < L; k++) {
					average[n] += data[n*L + k];
					average2 += data[n*L + k]*data[n*L + k];
				}
				average[n] /= (double)L;
				average2 /= (double)L;
				sigma[n] = sqrt(average2 - average[n]*average[n]);
			}
		}
	}	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockAverageVariance(vector<double>& average, vector<double>& variance) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockAverageVariance(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		average.clear();
		variance.clear();
		if (L == 1) {
			cout << "Warning [Stat :: BlockAverageVariance(vector<double>&, vector<double>&)]: the standard deviation is automatically set to zero";
			cout << " since each block has length 1." << endl;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average.push_back(data[n]);
				variance.push_back(0);
			}
		}
		else {
			double average2;
			for (unsigned int n = 0; n < Nblocks; n++) {
				average.push_back(0);
				average2 = 0;
				variance.push_back(0);
				
				for (unsigned int k = 0; k < L; k++) {
					average[n] += data[n*L + k];
					average2 += data[n*L + k]*data[n*L + k];
				}
				average[n] /= (double)L;
				average2 /= (double)L;
				variance[n] = average2 - average[n]*average[n];
			}
		}
	}	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Standard deviation (using expected mean), calculated for each block of data
double Stat :: BlockSigmaExp(vector<double>& sigma_exp) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockSigmaExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		sigma_exp.clear();
		if (L == 1) { // Blocks of length 1
			for (unsigned int n = 0; n < Nblocks; n++) sigma_exp[n] = abs(data[n]-exp_mean);
		}
		else  {
			for (unsigned int n = 0; n < Nblocks; n++) {
				sigma_exp.push_back(0);
				for (unsigned int k = 0; k < L; k++) {
					sigma_exp[n] += (data[n*L + k]-exp_mean)*(data[n*L + k]-exp_mean); // ## //
				}
				sigma_exp[n] = sqrt(sigma_exp[n]/(double)L);
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Quadratic deviation from expected mean, calculated for each block of data
double Stat :: BlockVarianceExp(vector<double>& variance_exp) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockVarianceExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		variance_exp.clear();
		if (L == 1) { // Blocks of length 1
			for (unsigned int n = 0; n < Nblocks; n++) variance_exp[n] = (data[n]-exp_mean)*(data[n]-exp_mean);
		}
		else {
			for (unsigned int n = 0; n < Nblocks; n++) {
				variance_exp.push_back(0);
				for (unsigned int k = 0; k < L; k++) {
					variance_exp[n] += (data[n*L + k]-exp_mean)*(data[n*L + k]-exp_mean);
				}
				variance_exp[n] /= (double)L;	
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


/// ****** ****** ****** ****** ******  Data Analysis - Progressive Methods (vector<double>) ****** ****** ****** ****** ****** ///


// Progressive <average> and <sigma of the average (sigma/sqrt(N-1)) of <data> vector (Equivalent to BlockProgAverage if Nblocks = 1)
double Stat :: ProgAverageSigma(vector<double>& prog_average, vector<double>& prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: ProgAverageSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		// Cumulative averages and errors
		prog_average.clear();
		prog_sigma.clear();
		vector<double> prog_average2;
		
		prog_average.push_back(data[0]);
		prog_average2.push_back(data[0]*data[0]);
		for (unsigned int i = 1; i < Nsample; i++) {
			prog_average.push_back(0);
			prog_average2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_average[i] += data[j];
				prog_average2[i] += data[j]*data[j];
			}
			prog_average[i] = prog_average[i] / (double)(i+1);
			prog_average2[i] = prog_average2[i] / (double)(i+1);
		}
		
		// Error: sigma with (N-1)
		prog_sigma.push_back(0);
		for (unsigned int i = 1; i < Nsample; i++) {
			prog_sigma.push_back( sqrt((prog_average2[i] - prog_average[i]*prog_average[i])/ (double)i) ); 
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Calculates the progressive average for each block (General form of ProgAverage(Nblocks = 1))
double Stat :: BlockProgAverageSigma(vector<double>& prog_average, vector<double>& prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgAverageSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_average, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_average.clear();
	prog_sigma.clear();
	
	block_average = 0;
	for (unsigned int k = 0; k < L; k++) {
		block_average += data[k];
	}
	block_average /= (double)L;
	total_sum = block_average;
	total_sum2 = block_average*block_average;
	
	prog_average.push_back(block_average); // Sum
	prog_sigma.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_average = 0;
		for (unsigned int k = 0; k < L; k++) {
			block_average += data[n*L + k];
		}
		block_average /= (double)L;
		total_sum += block_average;
		total_sum2 += block_average*block_average;
		prog_average.push_back(total_sum/(double)(n+1));
		prog_sigma.push_back( sqrt( (total_sum2/(double)(n+1)- prog_average[n]*prog_average[n])/ (double)n) ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgAverageVariance(vector<double>& prog_average, vector<double>& prog_variance) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgAverageVariance(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_average, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_average.clear();
	prog_variance.clear();
	
	block_average = 0;
	for (unsigned int k = 0; k < L; k++) {
		block_average += data[k];
	}
	block_average /= (double)L;
	total_sum = block_average;
	total_sum2 = block_average*block_average;
	
	prog_average.push_back(block_average); // Sum
	prog_variance.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_average = 0;
		for (unsigned int k = 0; k < L; k++) {
			block_average += data[n*L + k];
		}
		block_average /= (double)L;
		total_sum += block_average;
		total_sum2 += block_average*block_average;
		prog_average.push_back(total_sum/(double)(n+1));
		prog_variance.push_back( (total_sum2/(double)(n+1)- prog_average[n]*prog_average[n])/ (double)n ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgSigmaSigma(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigmaSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_average, block_average2, block_variance, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_sigma.clear();
	prog_error.clear();
	
	block_average = 0;
	block_average2 = 0;
	for (unsigned int k = 0; k < L; k++) {
		block_average += data[k];
		block_average2 += data[k]*data[k];
	}
	block_average /= (double)L;
	block_average2 /= (double)L;
	block_variance = block_average2 - block_average*block_average;
	total_sum2 = block_variance; // Sum of the blocks' variance
	total_sum = sqrt( block_variance ); // Sum of the blocks' standard deviations
	
	prog_sigma.push_back(total_sum); // Sum
	prog_error.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_average = 0;
		block_average2 = 0;
		for (unsigned int k = 0; k < L; k++) {
			block_average += data[n*L+k];
			block_average2 += data[n*L+k]*data[n*L+k];
		}
		block_average /= (double)L;
		block_average2 /= (double)L;
		block_variance = block_average2 - block_average*block_average;
		total_sum2 += block_variance; // Sum of the blocks' variance
		total_sum += sqrt( block_variance ); // Sum of the blocks' standard deviations
		
		prog_sigma.push_back(total_sum/(double)(n+1));
		prog_error.push_back( sqrt( (total_sum2/(double)(n+1)- prog_sigma[n]*prog_sigma[n])/ (double)n) ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgVarianceSigma(vector<double>& prog_variance, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgVarianceSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_average, block_average2, block_variance, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_variance.clear();
	prog_error.clear();
	
	block_average = 0;
	block_average2 = 0;
	for (unsigned int k = 0; k < L; k++) {
		block_average += data[k];
		block_average2 += data[k]*data[k];
	}
	block_average /= (double)L;
	block_average2 /= (double)L;
	block_variance = block_average2 - block_average*block_average;
	total_sum = block_variance; // Sum of the blocks' standard deviations
	total_sum2 = block_variance*block_variance; // Sum of the blocks' variance
	
	
	prog_variance.push_back(total_sum); // Sum
	prog_error.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_average = 0;
		block_average2 = 0;
		for (unsigned int k = 0; k < L; k++) {
			block_average += data[n*L+k];
			block_average2 += data[n*L+k]*data[n*L+k];
		}
		block_average /= (double)L;
		block_average2 /= (double)L;
		block_variance = block_average2 - block_average*block_average;
		total_sum += block_variance; // Sum of the blocks' standard deviations
		total_sum2 += block_variance*block_variance; // Sum of the blocks' variance
		
		prog_variance.push_back(total_sum/(double)(n+1));
		prog_error.push_back( sqrt( (total_sum2/(double)(n+1)- prog_variance[n]*prog_variance[n])/ (double)n) ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgSigmaSigmaExp(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigmaSigmaExp(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_variance_exp = 0, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_sigma.clear();
	prog_error.clear();

	for (unsigned int k = 0; k < L; k++) block_variance_exp += (data[k]-exp_mean)*(data[k]-exp_mean);
	block_variance_exp /= (double)L;
	total_sum2 = block_variance_exp; // Sum of the blocks' variance
	total_sum = sqrt( block_variance_exp ); // Sum of the blocks' standard deviations
	
	prog_sigma.push_back(total_sum); // Sum
	prog_error.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_variance_exp = 0;
		for (unsigned int k = 0; k < L; k++) block_variance_exp += (data[n*L+k]-exp_mean)*(data[n*L+k]-exp_mean);
		block_variance_exp /= (double)L;
		total_sum2 += block_variance_exp; // Sum of the blocks' variance
		total_sum += sqrt( block_variance_exp ); // Sum of the blocks' standard deviations
		
		prog_sigma.push_back(total_sum/(double)(n+1));
		prog_error.push_back( sqrt( (total_sum2/(double)(n+1)- prog_sigma[n]*prog_sigma[n])/ (double)n) ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgVarianceSigmaExp(vector<double>& prog_variance, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgVarianceSigmaExp(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	double block_variance_exp = 0, total_sum, total_sum2;
	// Cumulative averages and errors
	prog_variance.clear();
	prog_error.clear();

	for (unsigned int k = 0; k < L; k++) block_variance_exp += (data[k]-exp_mean)*(data[k]-exp_mean);
	block_variance_exp /= (double)L;
	total_sum = block_variance_exp; // Sum of the blocks' standard deviations
	total_sum2 = block_variance_exp*block_variance_exp; // Sum of the blocks' variance
	
	
	prog_variance.push_back(total_sum); // Sum
	prog_error.push_back(0);
	
	for (unsigned int n = 1; n < Nblocks; n++) {
		block_variance_exp = 0;
		for (unsigned int k = 0; k < L; k++) block_variance_exp += (data[n*L+k]-exp_mean)*(data[n*L+k]-exp_mean);
		block_variance_exp /= (double)L;
		total_sum += block_variance_exp; // Sum of the blocks' standard deviations
		total_sum2 += block_variance_exp*block_variance_exp; // Sum of the blocks' variance
		
		prog_variance.push_back(total_sum/(double)(n+1));
		prog_error.push_back( sqrt( (total_sum2/(double)(n+1)- prog_variance[n]*prog_variance[n])/ (double)n) ); // n = (N-1)
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

/* OLD
double Stat :: BlockProgAverageSigma(vector<double>& prog_average, vector<double>& prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	
	vector<double> average;
	BlockAverage(average);
	// Cumulative averages and errors
	prog_average.clear();
	prog_sigma.clear();
	vector<double> prog_average2;
	
	prog_average.push_back(average[0]);
	prog_average2.push_back(average[0]*average[0]);
	for (unsigned int i = 1; i < Nblocks; i++) {
		prog_average.push_back(0);
		prog_average2.push_back(0);
		for (unsigned int j = 0; j <= i; j++) {
			prog_average[i] += average[j];
			prog_average2[i] += average[j]*average[j];
		}
		prog_average[i] /= (double)(i+1);
		prog_average2[i] /= (double)(i+1);
	}
	
	// Error: sigma with (N-1)
	prog_sigma.push_back(0);
	for (unsigned int i = 1; i < Nblocks; i++) {
		prog_sigma.push_back( sqrt( (prog_average2[i] - prog_average[i]*prog_average[i])/ (double)i) ); 
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
*/

/* OLD
double Stat :: BlockProgAverageVariance(vector<double>& prog_average, vector<double>& prog_variance) {
	auto t1 = chrono::high_resolution_clock::now();
	
	vector<double> average;
	BlockAverage(average);
	
	// Cumulative averages and errors
	prog_average.clear();
	prog_variance.clear();
	vector<double> prog_average2;
	
	prog_average.push_back(average[0]);
	prog_average2.push_back(average[0]*average[0]);
	for (unsigned int i = 1; i < Nblocks; i++) {
		prog_average.push_back(0);
		prog_average2.push_back(0);
		for (unsigned int j = 0; j <= i; j++) {
			prog_average[i] += average[j];
			prog_average2[i] += average[j]*average[j];
		}
		prog_average[i] /= (double)(i+1);
		prog_average2[i] /= (double)(i+1);
	}
	
	// Error: sigma with (N-1)
	prog_variance.push_back(0);
	for (unsigned int i = 1; i < Nblocks; i++) {
		prog_variance.push_back( (prog_average2[i] - prog_average[i]*prog_average[i])/ (double)i ) ; 
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
*/

/* OLD
// Progressive sigma and its error (using each block's average) for each block of data
double Stat :: BlockProgSigmaSigma(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigma(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		vector<double> block_sigma;
		BlockSigma(block_sigma);
		
		prog_sigma.clear();
		prog_error.clear();
		vector<double> prog_sigma2;
		
		prog_sigma.push_back(block_sigma[0]);
		prog_sigma2.push_back(block_sigma[0]*block_sigma[0]);
		for (unsigned int i = 1; i < Nblocks; i++) {
			prog_sigma.push_back(0);
			prog_sigma2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_sigma[i] += block_sigma[j];
				prog_sigma2[i] += block_sigma[j]*block_sigma[j];
			}
			prog_sigma[i] /= (double)(i+1);
			prog_sigma2[i] /= (double)(i+1);
		}
	
		// Progressive sigma of the mean sigma ( = sigma_(measure)/sqrt(N-1) )
		prog_error.push_back(0); // First block's sigma is not computable, so it is set to zero
		for (unsigned int i = 1; i < Nblocks; i++) {
			prog_error.push_back( sqrt( (prog_sigma2[i] - prog_sigma[i]*prog_sigma[i])/ (double)i ) ); 
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
*/

/* OLD
// Progressive sigma and its error (considering an expected mean), calculated for each block of data
double Stat :: BlockProgSigmaSigmaExp(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigmaSigmaExp(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		vector<double> block_sigma_exp;
		BlockSigmaExp(block_sigma_exp);
	
		prog_sigma.clear();
		prog_error.clear();
		vector<double> prog_sigma2;
		
		prog_sigma.push_back(block_sigma_exp[0]);
		prog_sigma2.push_back(block_sigma_exp[0]*block_sigma_exp[0]);
		for (unsigned int i = 1; i < Nblocks; i++) {
			prog_sigma.push_back(0);
			prog_sigma2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_sigma[i] += block_sigma_exp[j];
				prog_sigma2[i] += block_sigma_exp[j]*block_sigma_exp[j];
			}
			prog_sigma[i] /= (double)(i+1);
			prog_sigma2[i] /= (double)(i+1);
		}
	
		// Progressive sigma of the mean sigma ( = sigma_(measure)/sqrt(N-1) )
		prog_error.push_back(0);
		for (unsigned int i = 1; i < Nblocks; i++) {
			prog_error.push_back( sqrt((prog_sigma2[i] - prog_sigma[i]*prog_sigma[i])/ (double)i) ); 
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
*/

/// ****** ****** ****** ****** ******  Data Analysis - Progressive Methods (last double&) ****** ****** ****** ****** ****** ///


double Stat :: BlockProgAverageSigma(double& last_prog_average, double& last_prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgAverage(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	if (Nblocks == 0) {
		cerr << "Error [Stat :: BlockProgAverage(double&, double&)]: Nblocks <= 1" << endl;
		exit(-1);
	}	
	else {
		last_prog_average = 0;
		double block_average, last_prog_average2 = 0;
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		for (unsigned int n = 0; n < Nblocks; n++) {
			block_average = 0;
			if (L == 1) block_average = data[n]; // Blocks of length 1
			else {
				for (unsigned int k = 0; k < L; k++) block_average += data[n*L + k];
				block_average /= (double)L;
			}
			last_prog_average += block_average; // Average of the block means
			last_prog_average2 += block_average*block_average;
		}
		last_prog_average /= (double)Nblocks;
		last_prog_average2 /= (double)Nblocks;
		last_prog_sigma = sqrt( (last_prog_average2 - last_prog_average*last_prog_average) / (double)(Nblocks-1) );
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgAverageVariance(double& last_prog_average, double& last_prog_variance) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgAverageVariance(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	if (Nblocks <= 1) {
		cerr << "Error [Stat :: BlockProgAverageVariance(double&, double&)]: Nblocks <= 1" << endl;
		exit(-1);
	}	
	else {
		last_prog_average = 0;
		double block_average, last_prog_average2 = 0;
		unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
		for (unsigned int n = 0; n < Nblocks; n++) {
			block_average = 0;
			if (L == 1) // Blocks of length 1
			{ 
				block_average = data[n];
			}
			else 
			{
				for (unsigned int k = 0; k < L; k++) block_average += data[n*L + k];
				block_average /= (double)L;
			}
			last_prog_average += block_average;
			last_prog_average2 += block_average*block_average;
		}
		last_prog_average /= (double)Nblocks;
		last_prog_average2 /= (double)Nblocks;
		last_prog_variance = (last_prog_average2 - last_prog_average*last_prog_average) / (double)(Nblocks-1) ;
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgSigmaSigma(double& last_prog_sigma, double& last_prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigma(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	if (Nblocks <= 1) {
		cerr << "Error [Stat :: BlockProgSigma(double&, double&)]: Nblocks <= 1" << endl;
		exit(-1);
	}	
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	if (L == 1) {
		cout << "Warning [Stat :: BlockProgSigma(double&, double&)]: the standard deviation is automatically set to zero";
		cout << " since each block has length 1." << endl;
		last_prog_sigma = 0;
		last_prog_error = 0;
	}
	else {
		last_prog_sigma = 0;
		double last_prog_sigma2 = 0;
		double block_average, block_average2;
		for (unsigned int n = 0; n < Nblocks; n++) {
			block_average = 0;
			block_average2 = 0;
			for (unsigned int k = 0; k < L; k++) {
				block_average += data[n*L+k];
				block_average2 += data[n*L+k]*data[n*L+k];
			}
			
			block_average /= (double)L;
			block_average2 /= (double)L;
			last_prog_sigma += sqrt(block_average2 - block_average*block_average);
			last_prog_sigma2 += block_average2 - block_average*block_average;
		}
		last_prog_sigma /= (double)Nblocks;
		last_prog_sigma2 /= (double)Nblocks;
		last_prog_error = sqrt( (last_prog_sigma2 - last_prog_sigma*last_prog_sigma) / (double)(Nblocks-1) );
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgVarianceSigma(double& last_prog_variance, double& last_prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgVarianceSigma(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	if (Nblocks <= 1) {
		cerr << "Error [Stat :: BlockProgVarianceSigma(double&, double&)]: Nblocks <= 1" << endl;
		exit(-1);
	}	
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	if (L == 1) {
		cout << "Warning [Stat :: BlockProgVarianceSigma(double&, double&)]: the standard deviation is automatically set to zero";
		cout << " since each block has length 1." << endl;
		last_prog_variance = 0;
		last_prog_error = 0;
	}
	else {
		last_prog_variance = 0;
		double last_prog_variance2 = 0;
		double block_average, block_average2;
		for (unsigned int n = 0; n < Nblocks; n++) {
			block_average = 0;
			block_average2 = 0;
			for (unsigned int k = 0; k < L; k++) {
				block_average += data[n*L+k];
				block_average2 += data[n*L+k]*data[n*L+k];
			}
			
			block_average /= (double)L;
			block_average2 /= (double)L;
			last_prog_variance += block_average2 - block_average*block_average; // Sum of block Variances
			last_prog_variance2 += (block_average2-block_average*block_average)*(block_average2-block_average*block_average); // Sum of block Variances squared
		}
		last_prog_variance /= (double)Nblocks;
		last_prog_variance2 /= (double)Nblocks;
		last_prog_error = sqrt( (last_prog_variance2 - last_prog_variance*last_prog_variance) / (double)(Nblocks-1) );
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: BlockProgSigmaSigmaExp(double& last_prog_sigma, double& last_prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (Nsample == 0) {
		cerr << "Error [Stat :: BlockProgSigma(double&, double&)]: Nsample = 0" << endl;
		exit(-1);
	}
	if (Nblocks <= 1) {
		cerr << "Error [Stat :: BlockProgSigma(double&, double&)]: Nblocks <= 1" << endl;
		exit(-1);
	}	
	unsigned int L = (unsigned int) Nsample/Nblocks; //Lenght of each block of data
	
	if (L == 1) {
		cout << "Warning [Stat :: BlockProgSigma(double&, double&)]: the standard deviation is automatically set to zero";
		cout << " since each block has length 1." << endl;
		last_prog_sigma = 0;
		last_prog_error = 0;
	}
	else {
		last_prog_sigma = 0;
		double block_sigma2, last_prog_sigma2 = 0;
		for (unsigned int n = 0; n < Nblocks; n++) {
			block_sigma2 = 0;
			for (unsigned int k = 0; k < L; k++) {
				block_sigma2 += (data[n*L+k]-exp_mean)*(data[n*L+k]-exp_mean);
			}
			block_sigma2 /= (double)L;
			last_prog_sigma += sqrt(block_sigma2);
			last_prog_sigma2 += block_sigma2;
		}
		last_prog_sigma /= (double)Nblocks;
		last_prog_sigma2 /= (double)Nblocks;
		last_prog_error = sqrt( (last_prog_sigma2 - last_prog_sigma*last_prog_sigma) / (double)(Nblocks-1) );
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

/// ****** ****** ****** ****** ****** ****** ******  Other Methods  ****** ****** ****** ****** ****** ****** ****** ///

double Stat :: Max() const {
	double max = data[0];
	for (unsigned int i = 1; i < data.size(); i++) {
		if (data[i] > max) max = data[i];
	}
	return max;
}

double Stat :: Max(unsigned int first, unsigned int last) const {
	if (first >= 0 && first < data.size() && last >= 0 && last < data.size()) {
		double max = data[first];
		for (unsigned int i = first+1; i < last; i++) {
			if (data[i] > max) max = data[i];
		}
		return max;
	}
	else {
		cerr << "Error [Stat :: Max]: invalid indices" << endl;
		exit(-1);
	}
	return 0;
}

double Stat :: Min() const {
	double min = data[0];
	for (unsigned int i = 1; i < data.size(); i++) {
		if (data[i] < min) min = data[i];
	}
	return min;
}

double Stat :: Min(unsigned int first, unsigned int last) const {
	if (first >= 0 && first < data.size() && last >= 0 && last < data.size()) {
		double min = data[first];
		for (unsigned int i = first+1; i < last; i++) {
			if (data[i] < min) min = data[i];
		}
		return min;
	}
	else {
		cerr << "Error [Stat :: Max]: invalid indices" << endl;
		exit(-1);
	}
	return 0;
}

void Stat :: PrintData() const {
	cout << "Data vector: " << endl;
	for (unsigned int i = 0 ; i < data.size(); i ++) {
		cout << " #" << i+1 << " -> " << data[i] << endl;
	}
	cout << "# # # # # # # # #" << endl;
	return;
}

void Stat :: Print(bool printdata = true) {
	cout << "Stat parameters: " << endl;
	cout << "Sample size: " << Nsample << endl;
	cout << "Number of Nblocks: " << Nblocks << endl;
	if (printdata) PrintData();
	return;
}

void Stat :: ClearData() {
	data.clear();
	Nsample = 0;
	Nblocks = 1;
	return;
}

void Stat :: ClearAll() {
	ClearData();
	rnd->SaveSeed();
	rnd = nullptr;
	return;
}
