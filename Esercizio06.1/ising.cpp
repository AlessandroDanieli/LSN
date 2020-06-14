#include <iostream>
#include <fstream>
#include <chrono>
#include <ctype.h>
#include <cmath>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include "ising.h"
using namespace std;

// Input folder:  /            (same as the C++ code)
// Output folder: OutputData   (except for the output configuration)
// Observables: U, C, M, X

Ising :: Ising(): Nspin(0), Nobs(4), Spin(nullptr), T(0) {
	ObsMean = new double[Nobs];
	ObsSigma = new double[Nobs];
}
Ising :: Ising(Random* _Rnd): Nspin(0), Nobs(4), Rnd(_Rnd), Spin(nullptr), T(0) {
	ObsMean = new double[Nobs];
	ObsSigma = new double[Nobs];
}
Ising :: ~Ising() {}

void Ising :: SetNspin(unsigned int _Nspin) {
	if (Spin == nullptr) {
		Nspin = _Nspin;
		Spin = new double[Nspin];
	}
	else {
		double* new_Spin = new double[_Nspin];
		for (unsigned int i = 0; i < (Nspin < _Nspin ? Nspin : _Nspin); i++) new_Spin[i] = Spin[i]; // Resize preserving the data
		delete[] Spin;
		Spin = new_Spin;
	}
	return;
}

void Ising :: SetT(double _T) {
	if (_T <= 0) {
		cerr << "Error [Ising :: SetT(double)]: invalid input T <= 0" << endl;
		exit(-1);
	}
	T = _T;
	beta = 1./T;
	return;
}

void Ising :: SetStat(unsigned int _Nblocks, unsigned int _Nsteps) {
	Nblocks = _Nblocks;
	Nsteps = _Nsteps;
	U_data.Set(Nblocks*Nsteps, Nblocks);
	S_data.Set(Nblocks*Nsteps, Nblocks);
}

// Set random spins
void Ising :: SetRandomSpin() {
	for (unsigned int i = 0; i < Nspin; i++) {
		if (Rnd->Uniform() >= 0.5) Spin[i] = 1;
		else Spin[i] = -1;
	}
	Measure();
	return;
}


// ****** ****** ****** ****** Flip spins, Equilibration (Automatic), Simulation ****** ****** ****** ****** //

double Ising :: Flip() { // Returns the Ratio of Acceptance
	unsigned int o; // Index of the spin to be flipped
	unsigned int Naccepted = 0;
	double delta_u;
	double q; // Ratio p(x_try)/p(x_prev), range of motion
	// Metropolis Algorithm
	if (AlgOption == 'm') { 
		for (unsigned int i = 0; i < Nspin; i++) {
			o = (unsigned int)(Rnd->Uniform()*Nspin);
			delta_u = 2.*Spin[o]*( J * (Spin[Pbc(o-1)] + Spin[Pbc(o+1)]) + Hext);
			q = exp(-beta*delta_u); // beta < +inf

			if (Rnd->Uniform() < q) {
				Spin[o] = -Spin[o];
				u += delta_u;
				s += 2.*Spin[o];
				Naccepted++;
			}
		}
		return Naccepted/(double)(Nspin); // Ratio of Acceptance
 	}
 	// Gibbs sampling
	else { 
		for (unsigned int i = 0; i < Nspin; i++) {
			o = (unsigned int)(Rnd->Uniform()*Nspin);
			delta_u = 2.*Spin[o]*( J*(Spin[Pbc(o-1)] + Spin[Pbc(o+1)]) + Hext); // E_f - E_i due to a possible o-th spin flip
			q = 1./(1.+exp(beta*delta_u));
			
			if (Rnd->Uniform() < q) {
				Spin[o] = -Spin[o];
				u += delta_u;
				s += 2.*Spin[o];
			}
		}
		return 1; // Ratio of Acceptance
	}
}

void Ising :: AutoEquilibration(unsigned int Ns, double tolerance) {
	cout << "  Automatic Equilibration: ";
	double av = 0, av2 = 0;
	unsigned int k = 0;
	do {
		for (unsigned int i = 0; i < Ns; i++) {
			Flip();
			av += u;
			av2 += u*u;	
		}
		av /= (double) Ns;
		av2 /= (double) Ns;
		k++;
	} while ( av2 - av*av > tolerance*tolerance && k < 2000); // Limiting the number of iterations to 2000
	cout << k << " attempts made" << endl;
	return;
}

/* Manual Equilibration
double Ising :: Equilibration() {
	auto t1 = chrono::high_resolution_clock::now();
	string option;
	cout << "*** Equilibration Cycle ***" << endl;
	cout << "> Start equilibration cycles (e or E) or proceed with the simulation (any key): ";
	cin >> option;
	if (option != "e" && option != "E") return 0;
	else {
		cout << "Insert the number of steps for each equilibration cycle: ";
		unsigned int Ncycles;
		cin >> Ncycles;
		
		do {
			for (unsigned int i = 0; i < Ncycles; i++) Flip();
			Measure();
			cout << endl << "> Total Energy U = " << u << endl;
			cout << "> Repeat the equilibration cycles (e or E) or proceed with the simulation (any key): ";
			cin >> option;
			if (option != "e" && option != "E") break;
			
		} while (true);
		cout << "*** End of Equilibration Cycle ***" << endl << endl;
		auto t2 = chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
	}
}
*/

double Ising :: Simulate() {
	cout << "  Simulation Started..." << endl;
	auto t1 = chrono::high_resolution_clock::now();
	
	for (unsigned int i = 0; i < Nsteps*Nblocks; i++) {
		Flip();
		AppendMeasure();
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}

// ****** ****** ****** ****** Measures: Output, computation ****** ****** ****** ****** //

void Ising :: Measure() {
	u = 0; // Internal Energy
	s = 0; // Total Spin
	for (unsigned int i = 0; i < Nspin; i++) {
		u -= J*Spin[i]*Spin[Pbc(i+1)] + Hext*Spin[i];
		s += Spin[i]; 
	}
	return;
}

void Ising :: AppendMeasure() {
	U_data.Append(u);
	S_data.Append(s);
	return;
}
void Ising :: ClearMeasures() {
	U_data.ClearData();
	S_data.ClearData();
	return;
}

void Ising :: ComputeU() {
	U_data.BlockProgAverageSigma(ObsMean[0], ObsSigma[0]);
	ObsMean[0] /= (double)Nspin;
	ObsSigma[0] /= (double)Nspin;
	return;
}
void Ising :: ComputeC() {
	U_data.BlockProgVarianceSigma(ObsMean[1], ObsSigma[1]); // <C> = <H^2>-<H>^2 with H total energy
	ObsMean[1] *= beta*beta/(double)Nspin;
	ObsSigma[1] *= beta*beta/(double)Nspin;
}
void Ising :: ComputeM() {
	S_data.BlockProgAverageSigma(ObsMean[2], ObsSigma[2]);
	ObsMean[2] /= (double)Nspin;
	ObsSigma[2] /= (double)Nspin;
	return;
}
void Ising :: ComputeX() {
	S_data.BlockProgVarianceSigma(ObsMean[3], ObsSigma[3]); // <X> = <S^2>-<S>^2 with S total spin
	ObsMean[3] *= beta/(double)Nspin;
	ObsSigma[3] *= beta/(double)Nspin;
}
void Ising :: ComputeObs(int obs) {
	if (obs == -1) { ComputeU(); ComputeC(); ComputeM(); ComputeX(); } // All 4 default observables U, C, M, X are computed
	else if (obs == 0) ComputeU();
	else if (obs == 1) ComputeC();
	else if (obs == 2) ComputeM();
	else if (obs == 3) ComputeX();
	else {
		cerr << "Error [Ising :: ComputeObs(int)]: invalid observable index" << endl;
		exit(-1);
	}
	return;
}

void Ising :: OutputMeasure(int obs) const {
	ofstream Out;
	if (obs == -1) { // All observables are computed and printed to the corresponding file
		for (unsigned int i = 0; i < Nobs; i++) {
			Out.open(("OutputData/obs"+to_string(i)+"_"+string(1,AlgOption)+".out").c_str(), ios::app);
			if (!Out.is_open()) {
				cerr << "Error [Ising :: OutputMeasure(int)]: cannot open the output file" << endl;
				exit(-1);
			}
			Out << Nspin << " " << T << " " << J << " " << Hext << " " << ObsMean[i] << " " << ObsSigma[i] << endl;
			Out.close();
		}
	}
	else if (obs >= 0 && abs(obs) < Nobs) { // Print single observable
		Out.open(("OutputData/obs"+to_string(obs)+"_"+string(1,AlgOption)+".out").c_str(), ios::app);
		if (!Out.is_open()) {
			cerr << "Error [Ising :: OutputMeasure(int)]: cannot open the output file" << endl;
			exit(-1);
		}
		Out << Nspin << " " << T << " " << J << " " << Hext << " " << ObsMean[obs] << " " << ObsSigma[obs] << endl;
		Out.close();
	}
	else {
		cerr << "Error [Ising :: OutputMeasure(int)]: invalid observable index" << endl;
		exit(-1);
	}
	return;
}

// ****** ****** ****** ****** Load/Output configurations ****** ****** ****** ****** //

bool Ising :: FileExists(const string& filename) const {
	struct stat buffer;   
	return (stat(filename.c_str(), &buffer) == 0); 
}

void Ising :: LoadConfig(const string filename) {
	if (FileExists(filename)) {
		cout << "*** Loading configuration from <" << filename << "> ***" << endl << endl;
		ifstream In(filename);
		if (In.is_open()) {
			string buff;
			In >> buff >> buff;
			unsigned int n = (unsigned int)stoi(buff);
			if (Spin == nullptr) {
				Nspin = n;
				Spin = new double[Nspin];
				cout << "A" << endl;
			}
			else if (n != Nspin) {
				Nspin = n;
				delete[] Spin;
				Spin = new double[Nspin];
				cout << "B" << endl;
			}
			// Now Nspin = Input dimension
			double s;
			for (unsigned int i = 0; i < Nspin; i++) {
				In >> s;
				Spin[i] = s;
			}
			In.close();
			Measure(); // Measure the actual values of U, S (total energy, total spin)
		}
		else {
			cerr << "Error [Ising :: LoadConfig(string)]: cannot open the file <" << filename << ">" << endl << endl;
			exit(-1);
		}
	}
	else {
		cerr << "Error [Ising :: LoadConfig(string)]: cannot open the input file <" << filename << ">" << endl << endl;
		exit(-1);
	}
	return;
}

void Ising :: OutputConfig(const string filename) const {
	cout << "*** Saving the final configuration on <" << filename << "> ***" << endl;
	ofstream Output(filename);
	if (Output.is_open()) {
		Output << "#Spin " << Nspin << endl;
		for (unsigned int i = 0; i < Nspin; i++) {
			Output << Spin[i] << endl;
		}
		Output.close();
	}
	else {
		cerr << "Error [Ising :: OutputConfig(const string)]: cannot open the output file" << endl;
		exit(-1);
	}
	return;
}

int Ising :: CleanOutput(int obs) const {
	int r = 0;
	if (obs == -1) {
		for (unsigned int i = 0; i < Nobs; i++) {
			r = system( ("rm -rf OutputData/obs"+to_string(i)+"_"+string(1,AlgOption)+".out").c_str() );
		}
	}
	else if (obs >= 0 && abs(obs) < Nobs) r = system( ("rm -rf OutputData/obs"+to_string(obs)+"_"+string(1,AlgOption)+".out").c_str() );
	else {
		cerr << "Error [Ising :: CleanOutput(int)]: invalid observable index" << endl;
		exit(-1);
	}
	return r;
}

unsigned int Ising :: Pbc(int i) { // Periodic Boundary Conditions
    if (i < 0) return (unsigned int) i + Nspin;
    else if (abs(i) >= Nspin) return (unsigned int) i - Nspin;
    else return i;
}

void Ising :: PrintParameters() const {
	cout << "*** System Parameters ***" << endl;
	cout << "> Number of spins = " << Nspin << endl;
	cout << "> Number of blocks = " << Nblocks << endl;
	cout << "> Number of steps of each block = " << Nsteps << endl;
	if (T >= 0) cout << "> Temperature = " << T << endl;
	cout << "> Exchange interaction = " << J << endl;
	cout << "> External field = " << Hext << endl;
	if (AlgOption == 'm') cout << "> Using Metropolis Algorithm" << endl;
	else if (AlgOption == 'g') cout << "> Using Gibbs Algorithm" << endl;
	else {
		cerr << "Error [Ising :: Input(const string)]: Invalid Algorithm option. Use 'M' ('m') to apply the Metropolis Algorithm or";
		cerr << " 'G' ('g') to apply Gibbs Algorithm. " << endl;
		exit(-1);
	}
	cout << endl;
	return;
}
