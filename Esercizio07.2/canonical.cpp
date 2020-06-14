#include "canonical.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
using namespace std;

Canonical :: Canonical(): d(3), Np(0), Nb(1), Ns(0), Nbins(0), Naccepted(0), Nattempted(0), X(nullptr), g(nullptr), g_final(nullptr) {}
Canonical :: Canonical(Random* r): d(3), Np(0), Nb(1), Ns(0), Nbins(0), Naccepted(0), Nattempted(0), X(nullptr), Rnd(r), g(nullptr), g_final(nullptr) {}

Canonical :: ~Canonical() {
	Rnd->SaveSeed();
	if (X != nullptr) {
		for (unsigned int i = 0; i < d; i++) delete[] X[i];
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
void Canonical :: SetNbins(unsigned int _Nbins) {
	ResizeBins(_Nbins);
	return;
}
void Canonical :: SetBlocksSteps(unsigned int _Nb, unsigned int _Ns) {
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
					P_old += (1./h - 0.5)/h; // ++++++ //
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
					P_new += (1./h - 0.5)/h; // ++++++ //
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
	
	// Reset g(r) bins
	for (unsigned int b = 0; b < Nbins; b++) g[b] = 0;

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
				g[(unsigned int) (sqrt(dr2)/Rcut*(double)Nbins)]++;
			}
		}
	}
	u *= 4.;
	p *= 16./Vol;
	p += Rho*T;
	
	// NORMALIZATION. Consider once each i-j interaction, then multiply by 2 the final values of g!
	double Norm;
	for (unsigned int b = 0; b < Nbins; b++) {
		Norm = (Rho*(double)Np*2./3.*M_PI* ( pow((double)(b+1)/(double)Nbins*Rcut, 3) - pow((double)(b)/(double)Nbins*Rcut, 3)));
		
		g[b] /= Norm;
	}
}

void Canonical :: SampleG() {
	double dr2, dX[d];
	
	for (unsigned int b = 0; b < Nbins; b++) g[b] = 0; // Reset
	
	for (unsigned int i = 0; i < Np-1; i++) { // N(N-1)/2 Interactions
		for (unsigned int j = i+1; j < Np; j++) {
			dr2 = 0;
			for (unsigned int k = 0; k < d; k++) {
				dX[k] = Pbc(X[i][k] - X[j][k]);
				dr2 += dX[k]*dX[k];
			}
			if (dr2 < Rcut*Rcut) {
				g[(unsigned int) ( sqrt(dr2)/Rcut*(double)Nbins)]++;
			}
		}
	}
	
	// NORMALIZATION. Consider once each i-j interaction, then multiply by 2 the final values of g!
	double Norm;
	for (unsigned int b = 0; b < Nbins; b++) {
		Norm = (Rho*(double)Np*2./3.*M_PI* ( pow((double)(b+1)/(double)Nbins*Rcut, 3) - pow((double)(b)/(double)Nbins*Rcut, 3)));
		
		g[b] /= Norm;
	}
}


double Canonical :: Simulate() {
	auto t1 = chrono::high_resolution_clock::now();
	cout << "> Simulation Started..." << endl;
	Nattempted = 0; 
	Naccepted = 0;
	
	double* gBlockSum = new double[Nbins];
	int r = system("rm -rf OutputData/output.gofr.0"); r*= 1;
	ofstream OutputMeanG("OutputData/output.gofr.0", ios::app);
	OutputMeanG << "Bins: " << Nbins << endl;
	OutputMeanG << "Nblocks: " << Nb << endl;
	OutputMeanG << "BoxEdge: " << L << endl;
	OutputMeanG << "CutoffRadius: " << Rcut << endl;
	
	vector<vector<double>> ProgG;
	for (unsigned int b = 0; b < Nbins; b++) {
		ProgG.emplace_back(vector<double>(Nb, 0));
	}
	
	
	for (unsigned int i = 0; i < Nb; i++) {
		cout << "[" << Ns*(i+1) << " steps]" << endl;
		
		if (i != 0) { // SAVE ON FILE THE BLOCK MEAN
			for (unsigned int b = 0; b < Nbins; b++) {
				gBlockSum[b] /= (double)Ns;
				ProgG[b][i-1] = gBlockSum[b];
			}
			for (unsigned int b = 0; b < Nbins; b++) {
				OutputMeanG << gBlockSum[b] << " ";
			}
			OutputMeanG << endl;
		}
		
		for (unsigned int b = 0; b < Nbins; b++) { // RESET THE BLOCK SUM
			gBlockSum[b] = 0;
		}
		
		for (unsigned int i = 0; i < Ns; i++) {
			Move();
			SampleG();
			for (unsigned int b = 0; b < Nbins; b++) { // Keep summing over the current block
				gBlockSum[b] += g[b];
			}
		}
	}
	
	// LAST BLOCK
	for (unsigned int b = 0; b < Nbins; b++) {
		gBlockSum[b] /= (double)Ns;
		ProgG[b][Nb-1] = gBlockSum[b];
	}
	for (unsigned int b = 0; b < Nbins; b++) {
		OutputMeanG << gBlockSum[b] << " ";
	}
	OutputMeanG.close();
	
	//*** OUTPUT: FINAL VALUES OF g
	double average, sigma;
	Stat G_Final;
	
	ofstream OutputFinalG("OutputData/output.gave.0");
	OutputFinalG << "Bins: " << Nbins << endl;
	OutputFinalG << "Nblocks: " << Nb << endl;
	for (unsigned int b = 0; b < Nbins; b++) {
		G_Final.Set(ProgG[b], Nb, Nb); // Set the b-th bin block averages, the sample size (number of blocks)
		G_Final.BlockProgAverageSigma(average, sigma);
		
		OutputFinalG << average << " " << sigma << endl;
	}
	OutputFinalG.close();
	
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

void Canonical :: ResizeBins(int _Nbins) {
	if (_Nbins <= 0) {
		cerr << "Error [Canonical :: ResizeBins(int)]: Invalid input Nbins <= 0" << endl;
		exit(-1);
	}
	else if (abs(_Nbins) != Nbins) {
		if (g == nullptr) {
			Nbins = _Nbins;
			g = new double[Nbins];
			g_final = new double[Nbins];
		}
		else {
			Nbins = _Nbins;
			delete[] g;
			delete[] g_final;
			g = new double[Nbins];
			g_final = new double[Nbins];
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

void Canonical :: OutputG(const string filename) const {
	ofstream Output(filename);
	if (Output.is_open()) {
		for (unsigned int b = 0; b < Nbins; b++) {
			Output << b+1 << " " << g[b] << endl;
		}
		Output.close();
	}
	else {
		cerr << "Error [Canonical :: OutputG(const string filename)]: cannot open the output file" << endl;
		exit(-1);
	}
	return;
}

double Canonical :: Pbc(double x) const {
	return x - L * rint(x/L);
}


void Canonical :: PrintParameters() const {
	cout << "*** System Parameters ***" << endl;
	cout << "> Number of particles = " << Np << endl;
	cout << "> Number of blocks = " << Nb << endl;
	cout << "> Number of steps of each block = " << Ns << endl;
	cout << "> Number of bins (g(r) sampling) = " << Nbins << endl;
	cout << "> Temperature = " << T << endl;
	cout << "> Volume = " << Vol << endl;
	cout << "> Box Edge = " << L << endl;
	cout << "> Cutoff Radius = " << Rcut << endl;
	cout << "> Range of motion (Metropolis uniform cubic sampling) = " << RangeOfMotion << endl;
	cout << endl;
	return;
}

void Canonical :: PrintG() const {
	unsigned int Sum = 0;
	for (unsigned int b = 0; b < Nbins; b++) {
		cout << "Bin[" << b+1 << "] -> " << g[b] << endl;
		Sum += g[b];
	}
	cout << "Total interactions (r < cutoff radius): " << Sum << endl;
	cout << endl;
}
