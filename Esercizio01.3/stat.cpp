#include <cmath>
#include "stat.h"
using namespace std;

Stat :: Stat() : rnd(nullptr), sample_size(0), blocks(1) {}

Stat :: Stat(Random* r, unsigned int _sample_size, unsigned int _blocks) {
	rnd = r;
	if (_sample_size < _blocks) {
		cerr << "Error [Stat :: Stat(Random* , int , int)]: sample_size < blocks" << endl;
		exit(-1);
	}
	else if (_sample_size % _blocks != 0) {
		cerr << "Error [Stat :: Stat(Random* , int , int)]: sample_size/blocks is not an integer" << endl;
		exit(-1);
	}
	else {
		sample_size = _sample_size;
		blocks = _blocks;
	}
}

Stat :: Stat(const vector<double>& x, unsigned int _sample_size, unsigned int _blocks): rnd(nullptr) {
	if (_sample_size < _blocks) {
		cerr << "Error [Stat :: Stat(const vector<double>& , int , int)]: sample_size < blocks" << endl;
		exit(-1);
	}
	else if (_sample_size % _blocks != 0) {
		cerr << "Error [Stat :: Stat(const vector<double>& , int , int)]: sample_size/blocks is not an integer" << endl;
		exit(-1);
	}
	else if (_sample_size > x.size()) {
		cerr << "Error [Stat :: Stat(const vector<double>& , int , int)]: sample_size exceeds data.size() " << endl;
		exit(-1);
	}
	else {
		sample_size = _sample_size;
		blocks = _blocks;
		for (unsigned int i = 0; i < x.size(); i++) {
			data.push_back(x[i]);
		}
	}
}

Stat :: Stat(const Stat& s) {
	rnd = s.GetRandomGen();
	data = s.GetData();
	sample_size = s.GetSampleSize();
	blocks = s.GetBlocks();
	exp_mean = s.GetExpMean();
}

Stat :: ~Stat() {
	data.clear();
	if (rnd != nullptr)	rnd->SaveSeed();
}

double Stat :: SetData(const vector<double>& x) { // Sets data and sample_size = data size, block lenght = 1
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	data = x;
	sample_size = x.size();
	blocks = x.size(); 
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: CopyData(const vector<double>& x) {
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	for (double el : x) data.push_back(el);
	sample_size = x.size();
	blocks = x.size();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Stat :: Set(const vector<double>& x, unsigned int _sample_size, unsigned int _blocks) {
	auto t1 = chrono::high_resolution_clock::now();
	data.clear();
	data = x;
	sample_size = _sample_size;
	blocks = _blocks;
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
	for (unsigned int i = 0; i < sample_size; i++) {
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
	for (unsigned int i = 0; i < sample_size; i++) {
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
	for (unsigned int i = 0; i < sample_size; i++) {
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
	for (unsigned int i = 0; i < sample_size; i++) {
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
	for (unsigned int i = 0; i < sample_size; i++) {
		data.push_back(rnd->Gauss(m, sigma));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Calculates the average of each block 
double Stat :: BlockAverage(vector<double>& average) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: BlockAverage(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}	
	else {
		unsigned int L = (unsigned int) sample_size/blocks; //Lenght of each block of data
		average.clear();
		for (unsigned int n = 0; n < blocks; n++) {
			average.push_back(0);
			if (L == 1) { // Blocks of length 1
				average[n] = data[n];
			}
			else {
				for (unsigned int k = 0; k < L; k++) {
					average[n] += data[n*L + k]; // ## //
				}
				average[n] /= L;
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Progressive average and sigma of <<data>> vector! -> (Equivalent to BlockProgAverage, setting blocks = 1)
double Stat :: ProgAverage(vector<double>& prog_average, vector<double>& prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: ProgAverage(vector<double>&, vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		// Cumulative averages and errors
		prog_average.clear();
		prog_sigma.clear();
		vector<double> prog_average2;
		
		prog_average.push_back(data[0]);
		prog_average2.push_back(data[0]*data[0]);
		for (unsigned int i = 1; i < sample_size; i++) {
			prog_average.push_back(0);
			prog_average2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_average[i] += data[j];
				prog_average2[i] += data[j]*data[j];
			}
			prog_average[i] = prog_average[i] / (i+1);
			prog_average2[i] = prog_average2[i] / (i+1);
		}
		
		// Error: sigma with (N-1)
		prog_sigma.push_back(0);
		for (unsigned int i = 1; i < blocks; i++) {
			prog_sigma.push_back( sqrt((prog_average2[i] - prog_average[i]*prog_average[i])/ i) ); 
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Calculates the progressive average for each block (General form of ProgAverage(blocks=1))
double Stat :: BlockProgAverage(vector<double>& prog_average, vector<double>& prog_sigma) {
	auto t1 = chrono::high_resolution_clock::now();
	
	vector<double> average;
	BlockAverage(average);
	
	// Cumulative averages and errors
	prog_average.clear();
	prog_sigma.clear();
	vector<double> prog_average2;
	
	prog_average.push_back(average[0]);
	prog_average2.push_back(average[0]*average[0]);
	for (unsigned int i = 1; i < blocks; i++) {
		prog_average.push_back(0);
		prog_average2.push_back(0);
		for (unsigned int j = 0; j <= i; j++) {
			prog_average[i] += average[j];
			prog_average2[i] += average[j]*average[j];
		}
		prog_average[i] /= (i+1);
		prog_average2[i] /= (i+1);
	}
	
	// Error: sigma with (N-1)
	prog_sigma.push_back(0);
	for (unsigned int i = 1; i < blocks; i++) {
		prog_sigma.push_back( sqrt( (prog_average2[i] - prog_average[i]*prog_average[i])/ i) ); 
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Quadratic deviation from expected mean, calculated for each block of data
double Stat :: BlockVarianceExp(vector<double>& variance_exp) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: BlockErrorExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		unsigned int L = (unsigned int) sample_size/blocks; //Lenght of each block of data
		variance_exp.clear();
		for (unsigned int n = 0; n < blocks; n++) {
			variance_exp.push_back(0);
			if (L == 1) { // Blocks of length 1
				variance_exp[n] = (data[n]-exp_mean)*(data[n]-exp_mean);
			}
			else {
				for (unsigned int k = 0; k < L; k++) {
					variance_exp[n] += (data[n*L + k]-exp_mean)*(data[n*L + k]-exp_mean); // ## //
				}
				variance_exp[n] /= L;
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Quadratic deviation from block average, calculated for each block of data
double Stat :: BlockVariance(vector<double>& error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: BlockErrorExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		vector<double> average;
		BlockAverage(average);
		unsigned int L = (unsigned int) sample_size/blocks; //Lenght of each block of data
		error.clear();
		for (unsigned int n = 0; n < blocks; n++) {
			error.push_back(0);
			if (L == 1) { // Blocks of length 1
				error[n] = (data[n]-average[n])*(data[n]-average[n]);
			}
			else {
				for (unsigned int k = 0; k < L; k++) {
					error[n] += (data[n*L + k]-average[n])*(data[n*L + k]-average[n]); // ## //
				}
				error[n] /= L;
			}
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Progressive sigma and its error (considering an expected variance), calculated for each block of data
double Stat :: BlockProgSigmaExp(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: BlockErrorExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		vector<double> variance_exp;
		BlockVarianceExp(variance_exp);
	
		prog_sigma.clear();
		prog_error.clear();
		vector<double> prog_sigma2;
		
		prog_sigma.push_back(variance_exp[0]);
		prog_sigma2.push_back(variance_exp[0]*variance_exp[0]);
		for (unsigned int i = 1; i < blocks; i++) {
			prog_sigma.push_back(0);
			prog_sigma2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_sigma[i] += variance_exp[j];
				prog_sigma2[i] += variance_exp[j]*variance_exp[j];
			}
			prog_sigma[i] /= (i+1);
			prog_sigma2[i] /= (i+1);
		}
	
		// Error: sigma with (N-1)
		prog_error.push_back(0);
		for (unsigned int i = 1; i < blocks; i++) {
			prog_error.push_back( sqrt((prog_sigma2[i] - prog_sigma[i]*prog_sigma[i])/ i) ); 
		}
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Progressive sigma and its error (considering each block's average), calculated for each block of data
double Stat :: BlockProgSigma(vector<double>& prog_sigma, vector<double>& prog_error) {
	auto t1 = chrono::high_resolution_clock::now();
	if (sample_size == 0) {
		cerr << "Error [Stat :: BlockErrorExp(vector<double>&)]: division by zero." << endl;
		exit(-1);
	}
	else {
		vector<double> variance;
		BlockVariance(variance);
		
		prog_sigma.clear();
		prog_error.clear();
		vector<double> prog_sigma2;
		
		prog_sigma.push_back(variance[0]);
		prog_sigma2.push_back(variance[0]*variance[0]);
		for (unsigned int i = 1; i < blocks; i++) {
			prog_sigma.push_back(0);
			prog_sigma2.push_back(0);
			for (unsigned int j = 0; j <= i; j++) {
				prog_sigma[i] += variance[j];
				prog_sigma2[i] += variance[j]*variance[j];
			}
			prog_sigma[i] /= (i+1);
			prog_sigma2[i] /= (i+1);
		}
	
		// Error: sigma with (N-1)
		prog_error.push_back(0);
		for (unsigned int i = 1; i < blocks; i++) {
			prog_error.push_back( sqrt( (prog_sigma2[i] - prog_sigma[i]*prog_sigma[i])/i ) ); 
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

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
		cout << " #" << i << " -> " << data[i] << endl;
	}
	cout << "# # # # # # # # #" << endl;
	return;
}

void Stat :: Print(bool printdata = true) {
	cout << "Stat parameters: " << endl;
	cout << "Sample size: " << sample_size << endl;
	cout << "Number of blocks: " << blocks << endl;
	if (printdata) PrintData();
	return;
}

void Stat :: ClearData() {
	data.clear();
	sample_size = 0;
	blocks = 1;
	return;
}

void Stat :: ClearAll() {
	ClearData();
	rnd->SaveSeed();
	rnd = nullptr;
	return;
}

