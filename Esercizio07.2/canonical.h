#ifndef __Canonical__
#define __Canonical__
#include <iostream>
#include <chrono>
#include "random.h"
#include "stat.h"
using namespace std;

class Canonical {
private:
	unsigned int d, Np, Nb, Ns, Nbins, Naccepted, Nattempted; // Dimensions, Number of particles, Number of blocks, Number of steps
	double **X; // Particle positions
	Random* Rnd;
	double T, beta, Vol, Rho, L, Rcut, RangeOfMotion; // Temperature, Thermodynamical beta, Volume, density, Box edge, cutoff radius, Metropolis range of motion (uniform cubic sampling)
	double u, p, *g, *g_final; // Instantaneous measures of energy and pressure, g(r) Radial distribution (2)
	
	
public:
	Canonical();
	Canonical(Random* Rnd);
	~Canonical();
	
	void ResizeNp(int _Np);
	void ResizeD(int _d);
	void ResizeBins(int _Nbins);
	
	// Set Methods
	void SetDim(unsigned int k);
	void SetNp(unsigned int k);
	void SetNbins(unsigned int k);
	void SetBlocksSteps(unsigned int _Nb, unsigned int _Ns);
	
	void SetT(double _T);
	void SetRho(double _Rho);
	void SetBox();
	void SetRcut(double _Rcut);
	void SetRangeOfMotion(double _RangeOfMotion);
	
	// Get Methods
	unsigned int GetDim() const { return d; }
	unsigned int GetNp() const { return Np; }
	unsigned int GetNb() const { return Nb; }
	unsigned int GetNs() const { return Ns; }
	unsigned int GetNbins() const { return Nbins; }
	double GetT() const { return T; }
	double GetBeta() const { return beta; }
	double GetRho() const { return Rho; }
	double GetVol() const { return Vol; }
	double GetL() const { return L; }
	double GetRcut() const { return Rcut; }
	double GetRatioAcc() const { return (double)Naccepted/(double)Nattempted; }
	double GetRangeOfMotion() const { return RangeOfMotion; }
	double GetInstU() const { return u; }
	double GetInstP() const { return p; }
	
	void AutoEquilibration(unsigned int Ns, double tolerance);
	double Equilibrium(unsigned int n_cycles, unsigned int move_steps);
	double Calibration(unsigned int Ncycles);
	
	void Move();
	double Simulate();
	void SampleG();
	void Measure();
	
	void LoadConfig(const string filename);
	void OutputConfig(const string filename) const;
	void OutputXYZ(unsigned int) const;
	//void OutputBlockMeasures(int obs) const;
	void OutputG(const string filename) const;
	double Pbc(double) const;
	
	void PrintParameters() const;
	void PrintG() const;
};

#endif
