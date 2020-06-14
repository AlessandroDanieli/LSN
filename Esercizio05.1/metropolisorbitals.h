#ifndef __Metropolisorbitals__
#define __Metropolisorbitals__
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "random.h"
#include "vector.h"
#include "orbital.h"
using namespace std;

class MetropolisOrbital : public Orbital {
private:
	Random* Rnd;
	Vector* Sample;
	double RatioAcc; // Acceptance ratio
	unsigned int d, Nsample; // Space dimension, Number of samples

public:
	MetropolisOrbital();
	MetropolisOrbital(Random* _Rnd, unsigned int _d, unsigned int _Nsample);
	MetropolisOrbital(Random* _Rnd, unsigned int _Nsample, Vector start);
	MetropolisOrbital(unsigned int _Z, unsigned int _n, unsigned int _l, int _m);
	MetropolisOrbital(Random* _Rnd, unsigned int _Nsample, Vector start, unsigned int _Z, unsigned int _n, unsigned int _l, int _m);
	~MetropolisOrbital();
	
	void SetRnd(Random* _Rnd) { Rnd = _Rnd; }
	void SetNsample(unsigned int _Nsample) { Nsample = _Nsample; }
	void SetStart(Vector& Start);
	
	Random* GetRnd() const { return Rnd; }
	unsigned int GetNsample() const { return Nsample; }
	unsigned int GetDim() const { return d; }
	double GetRatioAcc() const { return RatioAcc; }
	vector<Vector> GetSample() const;
	
	double Generate(unsigned int sampling, double& RangeOfMotion);
	double Calibrate(unsigned int sampling, double& RangeOfMotion, unsigned int CalibSteps);
	double GenerateUniform(double& RangeOfMotion); // Uniform Sampling
	double GenerateGaussian(double& RangeOfMotion); // Gaussian Sampling
	double CalibrateUniform(double& RangeOfMotion, unsigned int CalibSteps); // Uniform Sampling Calibration
	double CalibrateGaussian(double& RangeOfMotion, unsigned int CalibSteps); // Gaussian Sampling Calibration
	
	double MeanRadius() const;
	vector<double> GetRadius() const; // Returns the vector<double> containing all the radii of the Vectors generated
	void Print() const;
	void OutputSample(string filename, unsigned int first, unsigned int last, unsigned int step=1) const;
};

#endif
