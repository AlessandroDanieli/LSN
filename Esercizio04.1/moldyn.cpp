#include <cstdlib> 
#include <iostream>
#include <fstream>      
#include <cmath>
#include <chrono>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include "moldyn.h"
using namespace std;

// Output folder: OutputData
// Input configuration files: config.0, old.0, old.config
// Input parameters files: input.dat
// NOTE: the spatial dimensions are set by default to d = 3. In order to generalize this feature, part of the code should be modified.

MolDyn :: MolDyn(): X(nullptr), Xold(nullptr), V(nullptr), d(3), Np(0), Ns(0), Nprint(1), Nframe(1), Nmeasure(1) {}

MolDyn :: MolDyn(unsigned int _Ns, double _delta, unsigned int _Nprint, unsigned int _Nframe, unsigned int _Nmeasure): delta(_delta), d(3), Np(Np_default), Ns(_Ns), Nprint(_Nprint), Nframe(_Nframe), Nmeasure(_Nmeasure) {}

MolDyn :: ~MolDyn() {
	if (X != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] X[i];
		}
	}
	if (Xold != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] Xold[i];
		}
	}
	if (V != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] V[i];
		}
	}
}

void MolDyn :: SetParameters(unsigned int _Np, double _T, double _rho, double _rcut, double _dt, unsigned int _Ns, unsigned int _Nframe, unsigned int _Nmeasure, unsigned int _Nprint, unsigned int _Nblocks) {
	if (_Np < 2) {
		cerr << "Error [MolDyn :: SetParameters(double[...])]: Input number of particles < 2" << endl;
		exit(-1);
	}
	/*
	if (Np != _Np) {
		Clear(); // Deletes X, Xold, V in order to free the memory properly and initialize new dynamic arrays
		Np = _Np; // Input number of particles
	}
	*/
	InitializeArrays(_Np, 3);
	T = _T; // Input temperature
	rho = _rho; // Input volume
	rcut = _rcut; // Cutoff radius
	delta = _dt; // Time interval
	Ns = _Ns; // Number of steps
	Nframe = _Nframe; // Output frequency of frames
	Nmeasure = _Nmeasure; // Output frequency of measures
	Nblocks = _Nblocks;
	Nprint = _Nprint; // Current step info frequency
	vol = (double) Np/rho;
	L = pow(vol, 1.0/3.0);
	
	if (Nframe == 0) OutputF = false;
	if (Nmeasure == 0) OutputM = false;
	
	cout << "*** System parameters ***" << endl;
	cout << "> Number of particles = " << Np << endl;
	cout << "> Number of steps = " << Ns << endl ;
	cout << "> Density of particles = " << rho << endl;
	cout << "> Volume of the simulation box = " << vol << endl;
	cout << "> Edge of the simulation box = " << L << endl;
	cout << "> Time step = " << delta << endl;
	cout << "> Target Temperature = " << T << endl << endl;
	
	unsigned int Nm = (int)((double)(Ns)/(double)(Nmeasure));
	
	Prog_U.ClearData(); Prog_U.Set(Nm, Nblocks); 
	Prog_E.ClearData();	Prog_E.Set(Nm, Nblocks);
	Prog_K.ClearData();	Prog_K.Set(Nm, Nblocks); 
	Prog_T.ClearData();	Prog_T.Set(Nm, Nblocks); 

	return;
}


// Perform some cycles in order to reach an equilibrium temperature.
// n_cycles of <move_steps> steps each.
double MolDyn :: Equilibrium(unsigned int n_cycles, unsigned int move_steps) { // Scale speed in order to reach the target temperature T
	if (n_cycles*move_steps == 0) return 0;
	cout << "***  Equilibration cycle  ***" << endl;
	auto t1 = chrono::high_resolution_clock::now();
	//cout << "Velocity rescaling [Target Temperature = " << T << "]" << endl;
	
	double V_mid[Np][d], V2, ScaleFactor;
	ofstream EQ("OutputData/EQ.out");
	EQ << "Target: " << T << endl;
	
	for (unsigned int n = 0; n < n_cycles; n++) {
		Move(move_steps);
		MeasureT();
		cout << "  [Equilibration cycle " << n+1 << "] Measured temperature = " << E_T << endl;
		EQ << "[" << n+1 << "]	  " << E_T << endl;
		
		V2 = 0;
		for (unsigned int i = 1; i < Np; i++) {
			for (unsigned int k = 0; k < d; k++) {
				V_mid[i][k] = Pbc(X[i][k] - Xold[i][k])/delta;
				V2 += V_mid[i][k]*V_mid[i][k];
			}
		}
		
		V2 /= (double)(3*Np);
		ScaleFactor = sqrt(T/V2);
		for (unsigned int i = 1; i < Np; i++) {
			for (unsigned int k = 0; k < d; k++) {
				Xold[i][k] = Pbc(X[i][k] - V_mid[i][k]*delta*ScaleFactor);
			}
		}
	}
	EQ.close();
	cout << "> Temperature reached by the equilibration cycles: T = " << E_T << endl << endl;
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


// Optimized
void MolDyn :: Move(unsigned int n_cycles) { //Move particles with Verlet algorithm. V(t), r(t+dt)
	
	double **Xnew, **t;
	double F[d]; // Force
	double r[d], dr2, h;
	
	Xnew = new double*[Np];
	for (unsigned int i = 0; i < Np; i++) {
		Xnew[i] = new double[d];
	}
	for (unsigned int n = 0; n < n_cycles; n++) {
		for (unsigned int i = 0; i < Np; i++) {
			for (unsigned int k = 0; k < d; k++) {
				F[k] = 0;
			}
			for (unsigned int j = 0; j < Np; j++) {
				if (i != j) {
					dr2 = 0;
					for (unsigned int k = 0; k < d; k++) {
						r[k] = Pbc(X[i][k] - X[j][k]);
						dr2 += r[k]*r[k];
					}
					
					if (dr2 < rcut*rcut) { // Check whether the cutoff radius is greater than r_ij
						h = dr2*dr2*dr2;
						for (unsigned int k = 0; k < d; k++) {
							F[k] += r[k] / (h*dr2) * (2./h - 1.);
						}
					}
				}
			}
			for (unsigned int k = 0; k < d; k++) {
				Xnew[i][k] = Pbc(2*X[i][k] - Xold[i][k] + 24.*F[k]*delta*delta); // 24.0* added here and removed before from F[k] = ...
				V[i][k] = Pbc(Xnew[i][k] - Xold[i][k])/(2.0 * delta); // V(t) = (r(t+dt)-r(t-dt))/2dt
			}
		}
		t = Xold;
		Xold = X;
		X = Xnew;
		Xnew = t;
	}
	
	for (unsigned int i = 0; i < Np; i++) {
		delete[] Xnew[i];
	}
	delete[] Xnew;
	return;
}


double MolDyn :: Simulate() {
	cout << "***  Simulation  ***" << endl;
	auto t1 = chrono::high_resolution_clock::now();
	unsigned int frame_index = 1;
	ofstream OutU, OutK, OutE, OutT;
	/*
	if (Append == true) { // Append option
		OutU.open("OutputData/output_epot.dat", ios::app); 		OutK.open("OutputData/output_ekin.dat", ios::app);
		OutE.open("OutputData/output_etot.dat", ios::app);		OutT.open("OutputData/output_temp.dat", ios::app);
	}
	*/
	// Overwrite the existing files
	OutU.open("OutputData/output_epot.dat");	 OutK.open("OutputData/output_ekin.dat"); 
	OutE.open("OutputData/output_etot.dat");	 OutT.open("OutputData/output_temp.dat");
	
	if (!OutU.is_open() || !OutK.is_open() || !OutE.is_open() || !OutT.is_open()) {
		cerr << "[MolDyn :: Simulate]: An error occurred while opening the output files" << endl;
		exit(-1);
	}
	else {
		for (unsigned int s = 1; s <= Ns; s++) {
			Move();
			if (s % Nprint == 0) cout << "[" << s << " steps]"<< endl;
			if (OutputM && s % Nmeasure == 0) {
				Measure();
				OutputMeasure(OutU, OutK, OutE, OutT);
				AppendMeasure(); 	// Stat objects storing the instantaneous measured values
		 	}
		 	if (OutputF && s % Nframe == 0) {
		 		OutputFrame(frame_index);
		 		frame_index++;
		 	}
	  	}
		OutU.close(); OutK.close(); OutT.close(); OutE.close();
	}
	cout << endl;
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


void MolDyn :: Measure() { //Properties measurement
	double U = 0, K = 0, Uij, dr2, h;
	double r[d];

	// Potential Energy
	for (unsigned int i = 0; i < Np-1; i++) {
		for (unsigned int j = i+1; j < Np; j++) {
			dr2 = 0;
			for (unsigned int k = 0; k < d; k++) {
				r[k] = Pbc( Xold[i][k] - Xold[j][k] );
				dr2 += r[k]*r[k];
			}
			h = dr2*dr2*dr2;

			if (dr2 < rcut*rcut) {
				Uij = 4.0*(1./h - 1.)/h;
				U += Uij;
			}
		}
	}
	
	// Kinetic Energy
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			K += V[i][k]*V[i][k];
		}
	}
	K /= (double)Np;
	E_U = U/(double)Np; // Potential Energy per particle
	E_K = K*0.5; // Kinetic Energy per particle
	E_T = K/3.0; // Temperature
	E_E = (E_K+E_U); // Total Energy per particle
	return;
}

void MolDyn :: MeasureT() { // Measure T(t)
	E_T = 0;
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			E_T += V[i][k]*V[i][k];
		}
	}
	E_T /= (double)(3*Np);
	return;
}

void MolDyn :: AppendMeasure() {
	Prog_U.Append(E_U);
	Prog_K.Append(E_K);
	Prog_E.Append(E_E);
	Prog_T.Append(E_T);
}

void MolDyn :: OutputProgMeasure() { // Progressive averages of the observables
	ofstream Out;
	vector<double> ProgAverage, ProgSigma;
	
	Prog_U.BlockProgAverageSigma(ProgAverage, ProgSigma);
	Out.open("OutputData/output_epot_prog.dat"); 
	for (unsigned int i = 0; i < Nblocks; i++) {
		Out << ProgAverage[i] << " " << ProgSigma[i] << endl;
	}
	Out.close();
	Prog_K.BlockProgAverageSigma(ProgAverage, ProgSigma);
	Out.open("OutputData/output_ekin_prog.dat"); 
	for (unsigned int i = 0; i < Nblocks; i++) {
		Out << ProgAverage[i] << " " << ProgSigma[i] << endl;
	}
	Out.close();
	Prog_E.BlockProgAverageSigma(ProgAverage, ProgSigma);
	Out.open("OutputData/output_etot_prog.dat"); 
	for (unsigned int i = 0; i < Nblocks; i++) {
		Out << ProgAverage[i] << " " << ProgSigma[i] << endl;
	}
	Out.close();
	Prog_T.BlockProgAverageSigma(ProgAverage, ProgSigma);
	Out.open("OutputData/output_temp_prog.dat"); 
	for (unsigned int i = 0; i < Nblocks; i++) {
		Out << ProgAverage[i] << " " << ProgSigma[i] << endl;
	}
	Out.close();
	
	return;
}

void MolDyn :: OutputMeasure(ofstream& OutU, ofstream& OutK, ofstream& OutE, ofstream& OutT) {
	OutU << E_U  << endl;
	OutK << E_K  << endl;
	OutT << E_T << endl;
	OutE << E_E << endl;
	return;
}

void MolDyn :: OutputFrame(unsigned int index) {
	ofstream OutFrame;
	OutFrame.open("frames/config_" + to_string(index) + ".xyz");
	if (OutFrame.is_open()) {
		OutFrame << Np << endl << endl;
		for (unsigned int i = 0; i < Np; i++) {
			OutFrame << "LJ  ";
			for (unsigned int k = 0; k < d-1; k++) {
				OutFrame << Pbc(X[i][k]) << "   ";
			}
			OutFrame << Pbc(X[i][d-1]) << endl;
		}
		OutFrame.close();
	}
	else {
		cerr << "[MolDyn :: OutputFrame]: An error occurred while opening the " << index << "-th .xyz output file" << endl;
		exit(-1);
	}
	return;
}

double MolDyn :: OutputConfig(string OutNew, string OutOld) {
	auto t1 = chrono::high_resolution_clock::now();
	OutputConfig(OutNew);
	OutputOldConfig(OutOld);
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double MolDyn :: OutputConfig(string Out) { // Write final configuration
	auto t1 = chrono::high_resolution_clock::now();
	ofstream OutConf;
	OutConf.open(Out);
	if (OutConf.is_open()) {
		for (unsigned int i = 0; i < Np; i++) {
			for (unsigned int k = 0; k < d-1; k++) {
				OutConf << X[i][k]/L << "   ";
			}
			OutConf << X[i][d-1]/L << endl;
		}
		OutConf.close();
	}
	else {
		cerr << "[MolDyn :: OutputConf]: An error occurred while opening the output file" << endl;
		exit(-1);
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double MolDyn :: OutputOldConfig(string Out) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream OutConf;
	OutConf.open(Out);
	if (OutConf.is_open()) {
		for (unsigned int i = 0; i < Np; i++) {
			for (unsigned int k = 0; k < d-1; k++) {
				OutConf << Xold[i][k]/L << "   ";
			}
			OutConf << Xold[i][d-1]/L << endl;
		}
		OutConf.close();
	}
	else {
		cerr << "[MolDyn :: OutputOldConf]: An error occurred while opening the output file" << endl;
		exit(-1);
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}



double MolDyn :: SetRandomSpeed() {
	auto t1 = chrono::high_resolution_clock::now();
	double Vtot[d];
	
	if (X == nullptr) {
		cerr << "Error [MolDyn :: SetRandomSpeed()]: X is an empty array" << endl;
		exit(-1);
	}
	// Generates random initial speeds, such that Vtot = 0 and <K> = 3/2kT -> <v> = sqrt(3kT)
	/*
	if (V == nullptr) {
		V = new double*[Np];
		for (unsigned int i = 0; i < Np; i++) {
			V[i] = new double[d];
		}
	}
	*/
	for (unsigned int k = 0; k < d; k++) {
		V[0][k] = rand()/double(RAND_MAX) - 0.5; // Random values in (-0.5, 0.5)
		Vtot[k] = V[0][k];
	}
	for (unsigned int i = 1; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			V[i][k] = rand()/double(RAND_MAX) - 0.5;
			Vtot[k] += V[i][k];
		}
	}
	for (unsigned int k = 0; k < d; k++) {
		Vtot[k] /= (double) Np;
	}
	
	// Subtract the total velocity in order to eliminate the system drift (P_tot = 0 !)
	double V_2 = 0;
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			V[i][k] -= Vtot[k];
			V_2 += V[i][k]*V[i][k];
		}
	}
	V_2 /= (double) Np;
	double ScaleFactor = sqrt(3.0 * T / V_2);   // Velocity scale factor: Vtot^2 = 3kT
	
	// Old position generation (through random speed of each particle)
	////// Xold = new double*[Np];
	for (unsigned int i = 0; i < Np; i++) {
		//////  Xold[i] = new double[d];
		for (unsigned int k = 0; k < d; k++) {
			V[i][k] *= ScaleFactor;
			Xold[i][k] = Pbc(X[i][k] - V[i][k] * delta);
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

void MolDyn :: HalfStepSpeed(double** V_HS) {
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			V_HS[i][k] = Pbc(X[i][k] - Xold[i][k])/delta;
		}
	}
	return;
}

// Set the speeds computing one step, starting from two old configurations (old.0 and old.final)
double MolDyn :: SetSpeedInput() {
	auto t1 = chrono::high_resolution_clock::now();
	double Vtot[d];
	
	if (X == nullptr || Xold == nullptr) {
		cerr << "Error [MolDyn :: SetRandomSpeed()]: X or Xold are empty arrays" << endl;
		exit(-1);
	}
	/*
	if (V == nullptr) {
		V = new double*[Np];
		for (unsigned int i = 0; i < Np; i++) {
			V[i] = new double[d];
		}
	}
	*/
	Move(1);
	
	// Computation of the initial speeds at T=t+dt/2
	double** V_HS = new double*[Np];
	for (unsigned int i = 0; i < Np; i++) {
		V_HS[i] = new double[d];
	}
	HalfStepSpeed(V_HS);
	
	// Remove the drift (in this case, it should not be present, since the configurations are the outputs of the program. Computed as a precaution measure)
	for (unsigned int k = 0; k < d; k++) {
		Vtot[k] = V[0][k];
	}
	for (unsigned int i = 1; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			Vtot[k] += V[i][k];
		}
	}
	for (unsigned int k = 0; k < d; k++) {
		Vtot[k] /= (double) Np;
	}
	
	double V_2 = 0;
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			V[i][k] -= Vtot[k];
			V_2 += V[i][k]*V[i][k];
		}
	}
	V_2 /= (double) Np;
	double ScaleFactor = sqrt(3.0 * T / V_2);   // Velocity scale factor: Vtot^2 = 3kT
	
	// Old position generation through ScaleFactor correction
	for (unsigned int i = 0; i < Np; i++) {
		for (unsigned int k = 0; k < d; k++) {
			V[i][k] *= ScaleFactor;
			Xold[i][k] = Pbc(X[i][k] - V[i][k] * delta);
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


void MolDyn :: InitializeArrays(unsigned int New_Np, unsigned int New_d) {
	if (New_Np < 0) {
		cerr << "Error [MolDyn :: InitializeArrays(unsigned int, unsigned int)]: invalid new dimension for the arrays." << endl;
		exit(-1);
	}
	if (New_d < 0) {
		cerr << "Error [MolDyn :: InitializeArrays(unsigned int, unsigned int)]: invalid new space dimension." << endl;
		exit(-1);
	}
	if (New_Np == 0) {
		if (X == nullptr) {
			Np = 0; 
			return;
		}
		else {
			for (unsigned int i = 0; i < Np; i++) {
				delete[] X[i];
			}
			for (unsigned int i = 0; i < Np; i++) {
				delete[] Xold[i];
			}
			for (unsigned int i = 0; i < Np; i++) {
				delete[] V[i];
			}
			delete[] X;
			X = nullptr;
			delete[] Xold;
			Xold = nullptr;
			delete[] V;
			V = nullptr;
		}
	}
	else {
		if (X == nullptr) {
			Np = New_Np;
			d = New_d;
			X = new double*[Np];
			for (unsigned int i = 0; i < Np; ++i) {
				X[i] = new double[d];
			}
			Xold = new double*[Np];
			for (unsigned int i = 0; i < Np; ++i) {
				Xold[i] = new double[d];
			}
			V = new double*[Np];
			for (unsigned int i = 0; i < Np; ++i) {
				V[i] = new double[d];
			}
			return;
		}	
		else if (X != nullptr) {
			if (New_Np != Np) {
				for (unsigned int i = 0; i < Np; i++) {
					delete[] X[i];
				}
				for (unsigned int i = 0; i < Np; i++) {
					delete[] Xold[i];
				}
				for (unsigned int i = 0; i < Np; i++) {
					delete[] V[i];
				}
				delete[] X;
				delete[] Xold;
				delete[] V;
				
				Np = New_Np;
				d = New_d;
				
				X = new double*[Np];
				Xold = new double*[Np];
				V = new double*[Np];
				for (unsigned int i = 0; i < Np; ++i) {
					X[i] = new double[d];
				}
				for (unsigned int i = 0; i < Np; ++i) {
					Xold[i] = new double[d];
				}
				for (unsigned int i = 0; i < Np; ++i) {
					V[i] = new double[d];
				} 
				return;
			}
			else {
				if (New_d != d) {
					for (unsigned int i = 0; i < Np; i++) {
						delete[] X[i];
					}
					for (unsigned int i = 0; i < Np; i++) {
						delete[] Xold[i];
					}
					for (unsigned int i = 0; i < Np; i++) {
						delete[] V[i];
					}
					d = New_d;
					for (unsigned int i = 0; i < Np; ++i) {
						X[i] = new double[d];
					} 
					for (unsigned int i = 0; i < Np; ++i) {
						Xold[i] = new double[d];
					}
					for (unsigned int i = 0; i < Np; ++i) {
						V[i] = new double[d];
					} 
					return;
				}
				else return;
			}
		}
	}
}

double MolDyn :: SetConfig(string InputNew, string InputOld) {
	auto t1 = chrono::high_resolution_clock::now();
	SetConfig(InputNew); // Set X array reading <InputNew>
	SetOldConfig(InputOld); // Set Xold array reading <InputOld>
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}


// Used to check whether "old.0" and "old.final" exist
bool MolDyn :: FileExists(const string& name) {
	struct stat buffer;   
	return (stat(name.c_str(), &buffer) == 0); 
}

int MolDyn :: SetConfig(string filename) {
	bool InputExist = FileExists(filename);
	if (!InputExist) return -1;
	
	ifstream Input;
	Input.open(filename);
	if (!Input.is_open()) {
		cerr << "Error [MolDyn :: SetConfig]: invalid input filename" << endl;
		return -2;
	}

	for (unsigned int i = 0; i < Np; ++i) {
		for (unsigned int k = 0; k < d; k++) {
			Input >> X[i][k];
			X[i][k] *= L;
		}
	}
	Input.close();
	return 0;
}
int MolDyn :: SetOldConfig(string filename) {
	bool InputExist = FileExists(filename);
	if (!InputExist) return -1;
	
	ifstream Input;
	Input.open(filename);
	if (!Input.is_open()) {
		cerr << "Error [MolDyn :: SetOldConfig]: invalid input filename" << endl;
		return -2;
	}
	for (unsigned int i = 0; i < Np; ++i) {
		for (unsigned int k = 0; k < d; k++) {
			Input >> Xold[i][k];
			Xold[i][k] *= L;
		}
	}
	Input.close();
	return 0;
}


double MolDyn :: Pbc(double x) {
	return x - L*rint(x/L);
}

void MolDyn :: PrintX() const {
	for (unsigned int i = 0; i < Np; i++) {
		cout << "#" << i+1 << " > (";
		for (unsigned int k = 0; k < d-1; k++) {
			cout << X[i][k] << ", ";
		}
		cout << X[i][d-1] << ")" << endl;
	}
	cout << endl;
	return;
}
void MolDyn :: PrintXold() const {
	for (unsigned int i = 0; i < Np; i++) {
		cout << "#" << i+1 << " > (";
		for (unsigned int k = 0; k < d-1; k++) {
			cout << Xold[i][k] << ", ";
		}
		cout << Xold[i][d-1] << ")" << endl;
	}
	cout << endl;
	return;
}
void MolDyn :: PrintV() const {
	for (unsigned int i = 0; i < Np; i++) {
		cout << "#" << i+1 << " > (";
		for (unsigned int k = 0; k < d-1; k++) {
			cout << V[i][k] << ", ";
		}
		cout << V[i][d-1] << ")" << endl;
	}
	cout << endl;
	return;
}
void MolDyn :: PrintX(unsigned int first, unsigned int last) const {
	for (unsigned int i = first; i < last; i++) {
		cout << "#" << i+1 << " > (";
		for (unsigned int k = 0; k < d-1; k++) {
			cout << X[i][k] << ", ";
		}
		cout << X[i][d-1] << ")" << endl;
	}
	cout << endl;
}
void MolDyn :: PrintV(unsigned int first, unsigned int last) const {
	for (unsigned int i = first; i < last; i++) {
		cout << "#" << i+1 << " > (";
		for (unsigned int k = 0; k < d-1; k++) {
			cout << V[i][k] << ", ";
		}
		cout << V[i][d-1] << ")" << endl;
	}
	cout << endl;
	return;
}
void MolDyn :: PrintMeasures() const {
	cout << "Measures:" << endl;
	cout << "> Potential Energy = " << E_U << endl;
	cout << "> Kinetic Energy = " << E_K << endl;
	cout << "> Total Energy = " << E_E << endl;
	cout << "> Temperature = " << E_T << endl << endl;
	return;
}

void MolDyn :: Clear() {
	if (X != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] X[i];
		}
		delete[] X;
	}
	if (Xold != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] Xold[i];
		}
		delete[] Xold;
	}
	if (V != nullptr) {
		for (unsigned int i = 0; i < Np; i++) {
			delete[] V[i];
		}
		delete[] V;
	}
}
