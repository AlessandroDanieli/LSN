#ifndef __Stat__
#define __Stat__
#include <iostream>
#include <vector>
#include <cstring>
#include <chrono>
#include "random.h"
using namespace std;
// Block Average(Sigma/Variance): each block's average (and standard deviation/variance)
// Block Progressive Average(Sigma/Variance): mean of the entire sample (and standard deviation/variance of the mean = (<block_av^2>-<block_av>^2)/(N-1))

class Stat {
protected:
	Random* rnd;
	vector<double> data;
	unsigned int Nsample, Nblocks; // Number of Nblocks in which <data> is partitioned, sample size (not necessarily = data.size())
	double exp_mean, exp_sigma;
	
public:
	Stat();
	Stat(unsigned int _Nsample, unsigned int _Nblocks);
	Stat(Random* r, unsigned int _Nsample, unsigned int _Nblocks);
	Stat(const vector<double>& _data);
	Stat(const vector<double>& x, unsigned int _Nsample, unsigned int _Nblocks);
	Stat(const Stat& s);
	~Stat();
	double operator[](unsigned int i) const { return data[i]; }
	
	//## Set Methods ##//
	void SetRandomGen(Random* r) { rnd = r; }
	double SetData(const vector<double>&);
	double CopyData(vector<double> x);
	void SetSampleSize(unsigned int s) { Nsample = s; }
	void SetBlocks(unsigned int b) { Nblocks = b; }
	double Set(const vector<double>& x, unsigned int _Nsample, unsigned int _Nblocks);
	void Set(unsigned int _Nsample, unsigned int _Nblocks) { SetSampleSize(_Nsample); SetBlocks(_Nblocks); }
	void SetExpMean(double m) { exp_mean = m; }
	void SetExpSigma(double s) { exp_sigma = s; }
	void Append(double x);
	
	//## Get Methods ##//
	Random* GetRandomGen() const { return rnd; }
	vector<double> GetData() const { return data; }
	unsigned int GetSampleSize() const { return Nsample; }
	unsigned int GetDataSize() const { return data.size(); } // Might be different from size!
	unsigned int GetBlocks() const { return Nblocks; }
	double GetExpMean() const { return exp_mean; }
	double GetExpSigma() const { return exp_sigma; }
	
	//## Data Generation ##//
	double FillUniform();
	double FillUniform(double a, double b);
	double FillExp(double lambda);
	double FillLorentz(double m, double gamma);
	double FillGauss(double m, double sigma);
	
	//## Data Analysis ##//
	double Average(double& mean);
	double AverageSigma(double& mean, double& sigma);
	double AverageVariance(double& mean, double& variance);
	
	// *** Average with its Standard Deviation or Variance *** //
	double BlockAverage(vector<double>& average); 
	double BlockSigma(vector<double>& sigma);
	double BlockVariance(vector<double>& variance);
	double BlockAverageSigma(vector<double>& average, vector<double>& sigma);
	double BlockAverageVariance(vector<double>& average, vector<double>& sigma);
	// Progressive mean and sigma of the mean (sigma/sqrt(N-1)) with Nblocks = Nsample
	double ProgAverageSigma(vector<double>& prog_average, vector<double>& prog_sigma);
	// Progressive average and its error for each block of data
	double BlockProgAverageSigma(vector<double>& prog_average, vector<double>& prog_sigma);

	double BlockProgAverageVariance(vector<double>& prog_average, vector<double>& prog_variance);
	
	// Returning only the last progressive average / sigma / variance
	double BlockProgAverageSigma(double& last_prog_average, double& last_prog_sigma);
	double BlockProgAverageVariance(double& last_prog_average, double& last_prog_variance); 
	
	// ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** //
	// ***** Standard Deviation (or Variance) averages [w.r.t. block average or exp_mean] with its own Mean Standard Deviation  ***** //
	
	// Standard deviation (using the expected mean exp_mean) for each block of data
	double BlockSigmaExp(vector<double>& sigma_exp);
	// Variance (using the expected mean exp_mean) for each block of data
	double BlockVarianceExp(vector<double>& variance_exp);
	// Progressive sigma and its error (using each block's average) for each block of data
	double BlockProgSigmaSigma(vector<double>& prog_sigma, vector<double>& prog_error);
	// Progressive variance and its error (using each block's average) for each block of data
	double BlockProgVarianceSigma(vector<double>& prog_variance, vector<double>& prog_error);
	// Progressive sigma and its error (using the expected mean exp_mean) for each block of data
	double BlockProgSigmaSigmaExp(vector<double>& prog_sigma, vector<double>& prog_error);
	// Progressive variance and its error (using the expected mean exp_mean) for each block of data
	double BlockProgVarianceSigmaExp(vector<double>& prog_variance, vector<double>& prog_error);
	
	// Returning only the last progressive average / sigma / variance
	double BlockProgSigmaSigma(double& last_prog_sigma, double& last_prog_error);
	double BlockProgVarianceSigma(double& last_prog_variance, double& last_prog_error);
	double BlockProgSigmaSigmaExp(double& last_prog_sigma, double& last_prog_error);
	
	// ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** //
	
	double Max() const;
	double Min() const;
	double Max(unsigned int first, unsigned int last) const;
	double Min(unsigned int first, unsigned int last) const;
	
	void ClearData(); // data.clear()
	void ClearAll(); // clears all data (rnd = nullptr, ...)
	
	//## Output Methods ##//
	void PrintData() const;
	void Print(bool printdata);
};

#endif
