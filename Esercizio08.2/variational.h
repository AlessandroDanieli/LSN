#ifndef __Variational__
#define __Variational__
#include "random.h"
#include "stat.h"
using namespace std;

class Variational {
private:
	Stat H;
	Random* Rnd;
	double* Sample;
	double Mean, Sigma, RatioAcc, Start;
	double R_Mean, R_Sigma, R_Psi; // Ranges of motion for Mean, Sigma, Psi sampling (Metropolis algorithm - Uniform sampling)
	unsigned int Lblock, Nblocks, Nsample; // Total number of steps, number of blocks, number of points for each sampling of Psi^2
	double h; // <H> for the actual Sample
	bool Print = false;
	
	// Minimization parameters
	bool OptMinimize;
	unsigned int NsampleMinim; // Sample size (only for the minimize() function
	unsigned int TotalCycles, MoveMeanCycles, MoveSigmaCycles;
	double CycleScale, MoveMeanScale, MoveSigmaScale;
	
public:
	Variational();
	Variational(Random* _Rnd);
	
	bool FileExists(const string& filename);
	void Input(const string filename);
	void Resize(int _Nsample, int _Nblocks, int _Lblock);
	void ResizeSample(int _Nsample);
	void ResizeStat(int _Nblocks, int _Lblock);
	
	
	double Psi(double x, double m, double s) const;
	double Psi(double x) const;
	double d2_Psi(double x, double m, double s) const; // Second derivative (1D Laplacian)
	double d2_Psi(double x) const; // Second derivative (1D Laplacian)
	double H_Psi(double x, double m, double s) const;
	double H_Psi(double x) const;
	double p(double x, double m, double s) const;
	double p(double x) const;
	
	double V(double x) const;
	
	
	void SetStart(double x0);
	void SetMean(double m) { Mean = m; }
	void SetSigma(double s);
	void Set(int _Lblock, int _Nblocks, int _Nsample);
	void SetR_Psi(double r);
	void SetR_Mean(double r);
	void SetR_Sigma(double r);
	void SetPrint(bool b) { Print = b; }
	
	double GetExpH() { return h; }
	double GetMean() { return Mean; }
	double GetSigma() { return Sigma; }
	double GetRatioAcc() { return RatioAcc; }
	double* GetSample() { return Sample; }
	
	
	void GenerateUniform(double*& SampleNew, double m, double s);
	double GenerateUniform();
	
	double ComputeExpH(double* SampleNew, double m, double s);
	void ComputeExpH();
	
	void Calibrate();
	void Simulate();
	void Minimize();
	int MoveMean();
	void PrintMoveMean(string s, int C, double mleft, double mright, double hleft, double hright);
	double SafeSigma(double delta);
	int MoveSigma();
	void PrintMoveSigma(string s, int C, double sleft, double sright, double hleft, double hright);
	
	void OutputSample() const;
	void PrintSample() const;
	void PrintSample(double* SampleNew) const;
	
	void AppendH();
	void OutputProgH(vector<double> prog_average, vector<double> prog_sigma) const;
};

#endif
