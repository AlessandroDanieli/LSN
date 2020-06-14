#ifndef __Stat__
#define __Stat__
#include <iostream>
#include <vector>
#include <cstring>
#include <chrono>
#include "random.h"
#include "function.h"
using namespace std;

class Stat {
protected:
	Random* rnd;
	vector<double> data;
	unsigned int sample_size, blocks; // Number of blocks in which <data> is partitioned, sample size (not necessarily = data.size())
	double exp_mean;
	
public:
	Stat();
	Stat(Random* r, unsigned int _sample_size, unsigned int _blocks);
	Stat(const vector<double>& x);
	Stat(const vector<double>& x, unsigned int _sample_size, unsigned int _blocks);
	Stat(const Stat& s);
	~Stat();
	double operator[](unsigned int i) const { return data[i]; }
	//## Set Methods ##//
	void SetRandomGen(Random* r) { rnd = r; }
	double SetData(const vector<double>&);
	double CopyData(const vector<double>& x);
	void SetSampleSize(unsigned int s) { sample_size = s; }
	void SetBlocks(unsigned int b) { blocks = b; }
	double Set(const vector<double>&, unsigned int, unsigned int);
	void SetExpMean(double m) { exp_mean = m; }
	void Append(double x);
	//## Get Methods ##//
	Random* GetRandomGen() const { return rnd; }
	vector<double> GetData() const { return data; }
	unsigned int GetSampleSize() const { return sample_size; }
	unsigned int GetDataSize() const { return data.size(); } // Might be different from size!
	unsigned int GetBlocks() const { return blocks; }
	double GetExpMean() const { return exp_mean; }
	//## Data Generation and Manipulation Methods ##//
	double FillUniform();
	double FillUniform(double a, double b);
	double FillExp(double lambda);
	double FillLorentz(double m, double gamma);
	double FillGauss(double m, double sigma);
	
	double ProgAverage(vector<double>& prog_average, vector<double>& prog_sigma); // Average, Sigma of data vector (block length = 1)
	double BlockAverage(vector<double>& average);
	double BlockVariance(vector<double>& error);
	double BlockVarianceExp(vector<double>& variance_exp);
	double BlockProgAverage(vector<double>& prog_average, vector<double>& prog_sigma); 
	double BlockProgSigmaExp(vector<double>& prog_sigma, vector<double>& prog_error);
	double BlockProgSigma(vector<double>& prog_sigma, vector<double>& prog_error);
	
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

class StatIntegral : public Stat {
private:
	Function* f;
	vector<double> integral; // Vector composed of each block's integral value
	double xmin, xmax;
	
public:
	StatIntegral();
	StatIntegral(Random* r, unsigned int size, unsigned int bl, Function* g, double a, double b);
	~StatIntegral();
	
	//## Set Methods ##//
	void SetInterval(double a, double b) { xmin = a; xmax = b; }
	void SetFunction(Function*& g) { f = g; }
	
	//## Get Methods ##//
	Function* GetFunction() const { return f; }
	vector<double> GetIntegral() const { return integral; }	
	
	//## Other Methods ##//
	double FillLinear(double m, double q);
	double BlockIntegrate();
	double BlockIntegrate(Function* p);
	
	void PrintIntegral() const;
};

#endif
