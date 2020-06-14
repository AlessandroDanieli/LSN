#ifndef __Ising__
#define __Ising__
#include "random.h"
#include "stat.h"
#include <string>
using namespace std;

// Input folder:  /            (same as the C++ code)
// Output folder: OutputData   (except for the output configuration)
// Observables: U, C, M, X

class Ising {
private:
	unsigned int Nspin, Nblocks, Nsteps, Nobs;
	Random* Rnd;
	double* Spin;
	Stat U_data, S_data;
	double *ObsMean, *ObsSigma;
	
	char AlgOption; // Algorithm choice (Metropolis|Gibbs)
	double beta, T, J, Hext; // Parameters
	double u, s; // Instantaneous values of Internal Energy, Total Spin
	
public:
	Ising();
	Ising(Random* _Rnd);
	~Ising();
	
	// Set Methods
	void SetRandGen(Random* _Rnd) { Rnd = _Rnd; }
	void SetNspin(unsigned int _Nspin);
	void SetStat(unsigned int _Nblocks, unsigned int _Nsteps);
	void SetT(double _T);
	void SetJ(double _J) { J = _J; }
	void SetHext(double _H) { Hext = _H; }
	void SetAlgOption(char _opt) { AlgOption = _opt; }
	
	// Get Methods
	Random* GetRandGen() const { return Rnd; }
	unsigned int GetNspin() const { return Nspin; }
	unsigned int GetNblocks() const { return Nblocks; }
	unsigned int GetNsteps() const { return Nsteps; }
	double* GetSpin() const { return Spin; }
	double GetBeta() const { return beta; }
	double GetT() const { return T; }
	double GetJ() const { return J; }
	double GetHext() const { return Hext; }
	Stat GetUData() const { return U_data; }
	Stat GetSData() const { return S_data; }
	char GetAlgOption() const { return AlgOption; }
	
	
	// Initialization of the spin configuration
	void SetRandomSpin();
	// double Equilibration(); // Manual equilibration
	
	
	// Random flip of the spins
	double Flip();
	// Automatic equilibration over Ns cycles, using the variance of the energy and a threshold 
	void AutoEquilibration(unsigned int Ns, double tolerance);
	// Simulation using Metropolis or Gibbs sampling
	double Simulate(); // Flip(), AppendMeasure() repeated Nblocks*Nsteps times
	
	
	// Computation of observables
	void Measure();
	void AppendMeasure();
	void ClearMeasures();
	void ComputeU();
	void ComputeC();
	void ComputeM();
	void ComputeX();
	void ComputeObs(int obs = -1);	
	
	// Output and Input methods (Measures, configuration)
	bool FileExists(const string& filename) const;
	void OutputMeasure(int obs = -1) const;
	void LoadConfig(const string filename);
	void OutputConfig(const string filename) const;
	int CleanOutput(int obs = -1) const;
	
	unsigned int Pbc(int);
	void PrintParameters() const;
};

#endif

