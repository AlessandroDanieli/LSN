#include "canonical.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
using namespace std;

Canonical :: Canonical(): d(3), Np(0), Nb(1), Ns(0), Naccepted(0), Nattempted(0), X(nullptr) {}
Canonical :: Canonical(Random* r): d(3), Np(0), Nb(1), Ns(0), Naccepted(0), Nattempted(0), X(nullptr), Rnd(r) {}

Canonical :: ~Canonical() {
	Rnd->SaveSeed(); 
	if (X != nullptr) {
		for (unsigned int i = 0; i < d; i++) delete[] X[i];
		//delete[] X;
	}
}

void Canonical :: SetDim(unsigned int _d) {
	ResizeD(_d);
	return;
}
void Canonical :: SetNp(unsigned int _Np) {
	ResizeNp(_Np);
	return;
}

void Canonical :: SetBlocksSteps(unsigned int _Nb, unsigned int _Ns) {
	U.Set(_Nb*_Ns, _Nb);
	P.Set(_Nb*_Ns, _Nb);
	Nb = _Nb;
	Ns = _Ns;
	return;	
}

void Canonical :: SetT(double _T) {
	if (_T <= 0) {
		cerr << "Error [Canonical :: SetT(double)]: invalid input T <= 0" << endl;
		exit(-1);
	}
	T = _T;
	beta = 1./T;
	return;
}

void Canonical :: SetRho(double _Rho) {
	if (_Rho <= 0) {
		cerr << "Error [Canonical :: SetRho(double)]: invalid input Rho <= 0" << endl;
		exit(-1);
	}
	Rho = _Rho;
	return;
}

void Canonical :: SetBox() {
	Vol = (double)Np/Rho;
	L = pow(Vol, 1./3.);
	return;
}

void Canonical :: SetRcut(double _Rcut) {
	if (_Rcut <= 0) {
		cerr << "Error [Canonical :: SetRcut(double)]: invalid input Rcut <= 0" << endl;
		exit(-1);
	}
	Rcut = _Rcut;
	return;
}

void Canonical :: SetRangeOfMotion(double _RangeOfMotion) {
	if (_RangeOfMotion <= 0) {
		cerr << "Error [Canonical :: SetDelta(double)]: invalid input RangeOfMotion <= 0" << endl;
		exit(-1);
	}
	RangeOfMotion = _RangeOfMotion;
	return;
}

void Canonical :: AutoEquilibration(unsigned int Ns, double tolerance) {
	cout << "*** Equilibration ***" << endl;
	double av = 0, av2 = 0;
	unsigned int k = 0;
	do {
		av = 0;
		av2 = 0;
		for (unsigned int i = 0; i < Ns; i++) {
			Move();
			av += u;
			av2 += u*u;	
		}
		av /= (double) Ns;
		av2 /= (double) Ns;
		k++;
		cout << "  Attempt = " << k << "    U = " << av << "     Err% = " << sqrt(av2 - av*av)/abs(u)*100. << endl;
	} while ( ( k >= 6 && k < 20 && av2-av*av>tolerance*tolerance*(u*u) ) || k < 6); 
	cout << "  Automatic Equilibration: " << k << " attempts made" << endl;
	return;
}


void Canonical :: Move() {
	unsigned int o;
	
	//double p, energy_old, energy_new;
	double U_old, U_new, P_old, P_new, dr2, h, prob;
	double Xnew[d], dX[d];

	for (unsigned int i = 0; i < Np; i++) {
		o = (unsigned int)(Rnd->Uniform()*Np); // Randomly chosen particle index
		
		// ***** Old configuration and Energy ***** // // for (unsigned int k = 0; k < d; k++) Xold[k] = X[o][k];
		U_old = 0; // [o-th particle <-> all j!=o]
		P_old = 0; // [o-th particle <-> all j!=o] 
		for (unsigned int j = 0; j < Np; j++) {
			if (j != o) {
				dr2 = 0;
				for (unsigned int k = 0; k < d; k++) {
					dX[k] = Pbc(X[o][k] - X[j][k]);
					dr2 += dX[k]*dX[k];
				}
				
				if (dr2 < Rcut*Rcut) {
					h = pow(dr2, 3);
					U_old += (1./h - 1.)/h;
					P_old += (1./h - 0.5)/h;
				}
			}
		}
		U_old *= 4.;    // Old energy due to the interaction of [o] and [all j != o], for each [o] randomly chosen
		P_old *= 16./Vol;  
		P_old += Rho*T;
		
		// ***** Move the particle [o] ***** //
		for (unsigned int k = 0; k < d; k++) Xnew[k] = Pbc(X[o][k] + RangeOfMotion*(Rnd->Uniform() - 0.5));
		
		// ***** New configuration and Energy ***** //
		U_new = 0; P_new = 0; 
		for (unsigned int j = 0; j < Np; j++) {
			dr2 = 0;
			if (j != o) {
				for (unsigned int k = 0; k < d; k++) {
					dX[k] = Pbc(Xnew[k] - X[j][k]);
					dr2 += dX[k]*dX[k];
				}
				
				if (dr2 < Rcut*Rcut) { // Distance between o, j
					h = pow(dr2, 3);
					U_new += (1./h - 1.)/h;
					P_new += (1./h - 0.5)/h; 
				}
			}
		}
		U_new *= 4.; // New energy due to the interaction of [o] and [all j != o], for each [o] randomly chosen
		P_new *= 16./Vol;
		P_new += Rho*T;
		
		// ***** Metropolis probability test ***** //
		prob = exp(beta*(U_old-U_new)); 
		
		if (prob >= Rnd->Uniform()) {
			for (unsigned int k = 0; k < d; k++) X[o][k] = Xnew[k];
			u += (U_new - U_old);
			p += (P_new - P_old);
			Naccepted++;
		}
	}
	
	Nattempted += Np;
	return;
}

void Canonical :: Measure() {  // Unnecessary for u, p (already computed in Move())
	double dr2, h, dX[d];
	u = 0.0, p = 0.0;
	for (unsigned int i = 0; i < Np-1; i++) {
		for (unsigned int j = i+1; j < Np; j++) {
			// Distance (ij)
			dr2 = 0;
			for (unsigned int k = 0; k < d; k++) {
				dX[k] = Pbc(X[i][k] - X[j][k]);
				dr2 += dX[k]*dX[k];
			}

			if (dr2 < Rcut*Rcut) {
				h = pow(dr2, 3);
				u += (1./h - 1.)/h;
				p += (1./h - 0.5)/h;
			}
		}
	}
	u *= 4.;
	p *= 16./Vol;
	p += Rho*T;
	return;
}

double Canonical :: Simulate(const string filename) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream OutInst(filename);
	OutInst << "Cycles= " << Ns*Nb;
	
	cout << "> Simulation Started..." << endl;
	Nattempted = 0; 
	Naccepted = 0;
	for (unsigned int i = 0; i < Ns*Nb; i++) {
		if (i % Ns == 0) cout << "[" << i << " steps]" << endl;
		Move(); // Moves and updates U, P (Instantaneous)
		OutInst << endl << u/(double)Np << "   " << p;
		//SampleG();
		//AppendMeasure();
	}
	OutInst.close();
	cout << "Ratio of acceptance: " << (double)Naccepted/(double)Nattempted << endl;
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}

double Canonical :: Calibration(unsigned int Ncycles) {
	cout << "*** Calibration ***" << endl;
	Nattempted = 0; 
	Naccepted = 0;
	for (unsigned int i = 0; i < Ncycles; i++) {
		Move();
	}
	return (double)Naccepted/(double)Nattempted;
}


void Canonical :: LoadConfig(const string filename) {
	cout << "Loading configuration from <" << filename << ">..." << endl;
	ifstream Input(filename);
	unsigned int n;
	double x;
	if (Input.is_open()) {
		Input >> n;
		ResizeNp(n);
		for (unsigned int i = 0; i < Np; i++) for (unsigned int k = 0; k < d; k++) { Input >> x; X[i][k] = Pbc(x*L); }
		Input.close();
	}
	else {
		cerr << "Error [Canonical :: LoadConfig(const string)]: cannot open the input file" << endl;
		exit(-1);
	}
	Measure();
}


void Canonical :: AppendMeasure() {
	U.Append(u);
	P.Append(p);
	return;
}
void Canonical :: ClearMeasures() {
	U.ClearData();
	P.ClearData();
	return;
}

/*
void Canonical :: OutputMeasure(int obs) const {
	ofstream Out;
	if (obs == -1) { // All observables are computed and printed to the corresponding file
		for (unsigned int i = 0; i < Nobs; i++) {
			Out.open(("obs"+to_string(i)+".out").c_str(), ios::app);
			if (!Out.is_open()) {
				cerr << "Error [Canonical :: OutputMeasure(int)]: cannot open the output file" << endl;
				exit(-1);
			}
			Out << Np << " " << T << " " << Vol << " " << Rho << " " << Rcut << " " << ObsMean[i] << " " << ObsSigma[i] << endl;
			Out.close();
		}
	}
	else if (obs >= 0 && abs(obs) < Nobs) {
		Out.open(("obs"+to_string(obs)+".out").c_str(), ios::app);
		if (!Out.is_open()) {
			cerr << "Error [Canonical :: OutputMeasure(int)]: cannot open the output file" << endl;
			exit(-1);
		}
		Out << Np << " " << T << " " << Vol << " " << Rho << " " << Rcut << " " << ObsMean[obs] << " " << ObsSigma[obs] << endl;
		Out.close();
	}
	else {
		cerr << "Error [Canonical :: OutputMeasure(int)]: invalid observable index" << endl;
		exit(-1);
	}
	return;
}
*/
/*
int Canonical :: CleanOutput(int obs) const {
	int r = 0;
	if (obs == -1) {
		for (unsigned int i = 0; i < Nobs; i++) {
			r = system( ("rm -rf obs"+to_string(i)+".out").c_str() );
		}
	}
	else if (obs >= 0 && abs(obs) < Nobs) r = system( ("rm -rf obs"+to_string(obs)+".out").c_str() );
	else {
		cerr << "Error [Canonical :: CleanOutput(int)]: invalid observable index" << endl;
		exit(-1);
	}
	return r;
}
*/

void Canonical :: OutputBlockU() const {
	vector<double> mean, sigma;
	ofstream Out("OutputData/U_block.out");
	if (Out.is_open()) {
		U.BlockProgAverageSigma(mean, sigma);
		Out << "#Blocks " << Nb << endl;
		Out << "#BlockSteps " << Ns << endl;
		for (unsigned int n = 0; n < Nb; n++) {
			Out << mean[n] << " " << sigma[n] << endl;
		}
		cout << "OK" << endl;
		Out.close();
	}
	else {
		cerr << "Error [Canonical :: OutputBlockU()]: cannot open the output file" << endl;
		exit(-1);
	}
}

void Canonical :: OutputBlockP() const {
	vector<double> mean, sigma;
	ofstream Out("OutputData/P_block.out");
	if (Out.is_open()) {
		P.BlockProgAverageSigma(mean, sigma);
		Out << "#Blocks " << Nb << endl;
		Out << "#BlockSteps " << Ns << endl;
		for (unsigned int n = 0; n < Nb; n++) {
			Out << mean[n] << " " << sigma[n] << endl;
		}
		Out.close();
	}
	else {
		cerr << "Error [Canonical :: OutputBlockP()]: cannot open the output file" << endl;
		exit(-1);
	}
}

void Canonical :: ResizeNp(int Np_new) {
	if (Np_new < 0) {
		cerr << "Error [Canonical :: ResizeNp(int)]: Invalid input Np < 0" << endl;
		exit(-1);
	}
	else if (Np_new == 0) {
		if (X != nullptr) {
			for (unsigned int i = 0; i < Np; i++) delete[] X[i];
			delete[] X;
			X = nullptr;
			Np = 0;
		}
		else Np = 0;
	}
	else if (abs(Np_new) == Np) {
		if (X == nullptr) {
			X = new double*[Np];
			for (unsigned int i = 0; i < Np; i++) X[i] = new double[d];
		}
	}
	else {
		if (X == nullptr) {
			Np = Np_new;
			X = new double*[Np];
			for (unsigned int i = 0; i < Np; i++) X[i] = new double[d];
		}
		else {
			for (unsigned int i = 0; i < Np; i++) delete[] X[i];
			delete[] X;
			Np = Np_new;
			X = new double*[Np];
			for (unsigned int i = 0; i < Np; i++) X[i] = new double[d];
		}
	}
	return;
}

void Canonical :: ResizeD(int d_new) {
	if (d_new <= 0) {
		cerr << "Error [Canonical :: ResizeD(int)]: Invalid input d <= 0" << endl;
		exit(-1);
	}
	else if (abs(d_new) != d) {
		if (X == nullptr) {
			d = d_new;
			X = new double*[Np];
			for (unsigned int i = 0; i < Np; i++) X[i] = new double[d];
		}
		else {
			for (unsigned int i = 0; i < Np; i++) delete[] X[i];
			d = d_new;
			for (unsigned int i = 0; i < Np; i++) X[i] = new double[d];
		}
	}
	return;
}


// Reduced coordinates X_k / L
void Canonical :: OutputConfig(const string filename) const {
	ofstream Output;

	cout << "Outputting the particles' configuration to the file <" << filename << ">..." << endl << endl;
	Output.open("config.final");
	if (Output.is_open()) {
		Output << Np << endl;
		for (unsigned int i = 0; i < Np; i++) {
			for (unsigned int k = 0; k < d-1; k++) Output << X[i][k]/L << " ";
			Output << X[i][d-1]/L << endl;
		}
		Output.close();
	}
	else {
		cerr << "Error [Canonical :: OutputConfig(const string)]: cannot open the output file" << endl;
		exit(-1);
	}
	return;
}

// Real units (not reduced)
void Canonical :: OutputXYZ(unsigned int i_conf) const {
	ofstream Output("frames/config_" + to_string(i_conf) + ".xyz");
	if (Output.is_open()) {
		Output << Np << endl;
		for (unsigned int i = 0; i < Np; i++) {
			Output << "LJ ";
			for (unsigned int k = 0; k < d-1; k++) {
				Output << Pbc(X[i][k]) << " ";
			}
			Output << Pbc(X[i][d-1]) << endl;
		}
		Output.close();
	}
	else {
		cerr << "Error [Canonical :: ConfXYZ(unsigned int)]: cannot open the output file" << endl;
		exit(-1);
	}	
}

double Canonical :: Pbc(double x) const {
	return x - L * rint(x/L);
}


void Canonical :: PrintParameters() const {
	cout << "*** System Parameters ***" << endl;
	cout << "> Number of particles = " << Np << endl;
	cout << "> Number of blocks = " << Nb << endl;
	cout << "> Number of steps of each block = " << Ns << endl;
	cout << "> Temperature = " << T << endl;
	cout << "> Volume = " << Vol << endl;
	cout << "> Box Edge = " << L << endl;
	cout << "> Cutoff Radius = " << Rcut << endl;
	cout << "> Range of motion (Metropolis uniform cubic sampling) = " << RangeOfMotion << endl;
	cout << endl;
	return;
}
