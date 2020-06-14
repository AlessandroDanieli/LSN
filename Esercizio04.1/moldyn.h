#ifndef __MolDyn__
#define __MolDyn__
#include <vector>
#include "stat.h"
using namespace std;

class MolDyn {
private:
	// Position (actual and older), Velocity
	double **X = nullptr, **Xold = nullptr, **V = nullptr; 
	// Energy, temperature, volume, density, box edge length, cutoff radius, delta t
	double E, T, vol, rho, L, rcut, delta;
	// Estimate (per particle) of Potential energy, Kinetic energy, Total energy, and Temperature
	double E_U, E_K, E_E, E_T; 
	// Dimension (=3), Number of particles, Number of steps, Number of steps between each update, Number of steps between each frame output
	unsigned int d, Np, Ns, Nprint, Nframe, Nmeasure, Nblocks;
	// Default cubic lattice
	const unsigned int Np_default = 108;
	bool OutputF = true, OutputM = true;
	Stat Prog_U, Prog_K, Prog_E, Prog_T;
	
public:
	MolDyn();
	MolDyn(unsigned int _Ns, double _delta, unsigned int _Nprint, unsigned int _Nframe, unsigned int _Nmeasure);
	~MolDyn();
	
	//## Set Methods ##//
	void SetDim(unsigned int k) { d = k; }
	void SetNp(unsigned int k) { Np = k; }
	void SetNs(unsigned int k) { Ns = k; }
	void SetNprint(unsigned int k) { Nprint = k; }
	void SetNframe(unsigned int k) { Nframe = k; }
	void SetDelta(double k) { delta = k; }
	void SetOutputFrame(bool b) { OutputF = b; }
	void SetOutputMeasure(bool b) { OutputM = b; }
	void InitializeArrays(unsigned int New_Np, unsigned int New_d); // Initializes X, Xold, V (allocates/deallocates memory)
	//## Get Methods ##//
	unsigned int GetDim() const { return d; }
	unsigned int GetNp() const { return Np; }
	unsigned int GetNs() const { return Ns; }
	unsigned int GetNprint() const { return Nprint; }
	unsigned int GetNframe() const { return Nframe; }
	double GetDelta() const { return delta; }
	double GetT() const { return T; }
	double GetEstT() const { return E_T; }
	double GetEstU() const { return E_U; }
	double GetEstK() const { return E_K; }
	double GetEstE() const { return E_E; }
	bool GetOutputF() const { return OutputF; }
	
	void SetParameters(unsigned int _Np, double _T, double _rho, double _rcut, double _dt, unsigned int _Ns, unsigned int _Nfr, unsigned int _Nm, unsigned int _Npr, unsigned int _Nblocks);
	
	bool FileExists(const string& name);
	int SetConfig(string filename);
	int SetOldConfig(string filename);
	double SetConfig(string InputNew, string InputOld);

	double SetRandomSpeed(); // Prepare initial velocities, such that Vtot = 0 and <K> = 3/2kT -> <v> = sqrt(3kT)
	double SetSpeedInput();  // Prepare initial velocities, starting from two old configurations. Temperature correction is applied
		
	double Equilibrium(unsigned int n_cycles, unsigned int move_steps);
	void Move(unsigned int n_cycles = 1);
	void HalfStepSpeed(double** V_HS);
	//void Equilibrium(unsigned int n_cycles);
	
	void Measure();
	void MeasureT();
	double Simulate();
	double Pbc(double);
	
	void AppendMeasure();
	void OutputMeasure(ofstream& OutU, ofstream& OutK, ofstream& OutE, ofstream& OutT);
	void OutputProgMeasure();
	void OutputFrame(unsigned int index);
	double OutputConfig(string Out);
	double OutputOldConfig(string Out);
	double OutputConfig(string OutNew, string OutOld);
	void PrintX() const;
	void PrintXold() const;
	void PrintV() const;
	void PrintX(unsigned int first, unsigned int last) const;
	void PrintV(unsigned int first, unsigned int last) const;
	void PrintMeasures() const;
	void Clear();
};

#endif
