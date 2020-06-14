#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "variational.h"

Variational :: Variational(): Rnd(nullptr), Sample(nullptr), Mean(0), Sigma(1), Start(0), Lblock(0), Nblocks(0), Nsample(0) {}
Variational :: Variational(Random* _Rnd): Rnd(_Rnd), Sample(nullptr), Mean(0), Sigma(1), Start(0), Lblock(0), Nblocks(0), Nsample(0) {
}

double Variational :: Psi(double x, double m, double s) const {
	double diff = (x-m)/s;
	double sum = (x+m)/s;
	return exp(-diff*diff/2.) + exp(-sum*sum/2.);
}
double Variational :: Psi(double x) const {
	double diff = (x-Mean)/Sigma;
	double sum = (x+Mean)/Sigma;
	return exp(-diff*diff/2.) + exp(-sum*sum/2.);
}

double Variational :: d2_Psi(double x, double m, double s) const { // Laplacian (d/dx)^2 Psi(x)
	double diff = (x-m)/s;
	double sum = (x+m)/s;
	return ((diff*diff-1)*exp(-diff*diff/2.) + (sum*sum-1)*exp(-sum*sum/2.))/(s*s);
}
double Variational :: d2_Psi(double x) const { // Laplacian (d/dx)^2 Psi(x)
	double diff = (x-Mean)/Sigma;
	double sum = (x+Mean)/Sigma;
	return ((diff*diff-1)*exp(-diff*diff/2.) + (sum*sum-1)*exp(-sum*sum/2.))/(Sigma*Sigma);
}


double Variational :: H_Psi(double x, double m, double s) const {
	return -d2_Psi(x, m, s)/2. + V(x)*Psi(x, m, s);
}
double Variational :: H_Psi(double x) const {
	return -d2_Psi(x)/2. + V(x)*Psi(x);
}

double Variational :: p(double x, double m, double s) const { // p = |Psi|^2
	double t = Psi(x, m, s);
	return t*t;
}
double Variational :: p(double x) const { // p = |Psi|^2
	double t = Psi(x);
	return t*t;
}

void Variational :: ResizeSample(int _Nsample) {
	if (_Nsample <= 0) {
		cerr << "Error [Variational :: ResizeSample(int)]: Invalid input < 0" << endl;
		exit(-1);
	}
	else if (abs(_Nsample) == Nsample) {
		if (Sample == nullptr) {
			Sample = new double[Nsample];
		}
	}
	else {
		if (Sample == nullptr) {
			Nsample = _Nsample;
			Sample = new double[Nsample];
		}
		else {
			delete[] Sample;
			Nsample = _Nsample;
			Sample = new double[Nsample];
		}
	}
	return;
}

void Variational :: ResizeStat(int _Nblocks, int _Lblock) {
	if (_Lblock <= 0 || _Nblocks <= 0) {
		cerr << "Error [Variational :: ResizeStat(int, int)]: Invalid input < 0" << endl;
		exit(-1);
	}
	Lblock = _Lblock;
	Nblocks = _Nblocks;
	H.Set(Lblock*Nblocks, Nblocks); // Data size, Number of blocks
	return;
}

void Variational :: Resize(int _Nsample,  int _Nblocks, int _Lblock) {
	if (_Lblock <= 0 || _Nblocks <= 0 || _Nsample <= 0) {
		cerr << "Error [Variational :: Resize(int, int, int)]: Invalid input < 0" << endl;
		exit(-1);
	}
	else if (abs(_Nsample) == Nsample) {
		if (Sample == nullptr) {
			Sample = new double[Nsample];
		}
	}
	else {
		if (Sample == nullptr) {
			Nsample = _Nsample;
			Sample = new double[Nsample];
		}
		else {
			delete[] Sample;
			Nsample = _Nsample;
			Sample = new double[Nsample];
		}
	}
	Lblock = _Lblock;
	Nblocks = _Nblocks;
	H.Set(Lblock*Nblocks, Nblocks); // Data size, Number of blocks
	return;
}
void Variational :: Set(int _Lblock, int _Nblocks, int _Nsample) {
	Resize(_Lblock, _Nblocks, _Nsample);
	return;	
}
void Variational :: SetStart(double x0) {
	if (Sample == nullptr) {
		cerr << "Error [Variational :: SetStart(double)]: Sample array not initialized" << endl;
		exit(-1);
	}
	else {
		// Sample[0] = x0;
		Start = x0;
	}
	return;
}
void Variational :: SetSigma(double s) {
	if (s <= 0) {
		cerr << "Error [Variational :: SetSigma(double)]: Invalid input standard deviation" << endl;
		exit(-1);
	}
	Sigma = s;
	return;
}
void Variational :: SetR_Psi(double r) {
	if (r <= 0) {
		cerr << "Error [Variational :: SetR_Psi(double)]: Invalid input radius" << endl;
		exit(-1);
	}
	R_Psi = r;
	return;
}
void Variational :: SetR_Mean(double r) {
	if (r <= 0) {
		cerr << "Error [Variational :: SetR_Mean(double)]: Invalid input radius" << endl;
		exit(-1);
	}
	R_Mean = r;
	return;
}
void Variational :: SetR_Sigma(double r) {
	if (r <= 0) {
		cerr << "Error [Variational :: SetR_Sigma(double)]: Invalid input radius" << endl;
		exit(-1);
	}
	R_Sigma = r;
	return;
}

double Variational :: V(double x) const {
	int c = 1;
	if (c == 0) return x*x/2.;
	else if (c==1) {
		double t = x*x;
		return (t - 2.5)*t;
	}
	else return 0;
}



double Variational :: GenerateUniform() {
	auto t1 = chrono::high_resolution_clock::now();
	
	unsigned int Naccepted = 0;
	double q; // Ratio p(x_try)/p(x_prev), range of motion
	double x_try;
	
	Sample[0] = Start;
	for (unsigned int i = 1; i < Nsample; i++) {
		x_try = Sample[i-1] + Rnd->Uniform(-R_Psi, R_Psi);
		q = p(x_try)/p(Sample[i-1]);
		
		if (Rnd->Uniform() < q) {
			Sample[i] = x_try; 
			Naccepted++;
		}
		else {
			Sample[i] = Sample[i-1];
		}
	}
	RatioAcc = Naccepted/(double)(Nsample-1);
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}



void Variational :: GenerateUniform(double*& SampleNew, double m, double s) {
	//double* SampleNew = new double[Nsample];
	SampleNew[0] = Start;
	
	unsigned int Naccepted = 0;
	double q; // Ratio p(x_try)/p(x_prev), range of motion
	double x_try;
	
	for (unsigned int i = 1; i < Nsample; i++) {
		x_try = SampleNew[i-1] + Rnd->Uniform(-R_Psi, R_Psi);
		q = p(x_try, m, s)/p(SampleNew[i-1], m, s);
		
		if (Rnd->Uniform() < q) {
			SampleNew[i] = x_try; 
			Naccepted++;
		}
		else {
			SampleNew[i] = SampleNew[i-1];
		}
	}
	RatioAcc = Naccepted/(double)(Nsample-1);
	
	return;
}




void Variational :: ComputeExpH() {
	h = 0;
	for (unsigned int i = 0; i < Nsample; i++) {
		h += H_Psi(Sample[i])/Psi(Sample[i]);
	}
	h /= (double)Nsample;
}

double Variational :: ComputeExpH(double* SampleNew, double m, double s) {
	double h_new = 0;
	for (unsigned int i = 0; i < Nsample; i++) {
		h_new += H_Psi(SampleNew[i], m, s)/Psi(SampleNew[i], m, s);
	}
	h_new /= (double)Nsample;
	return h_new;
}

double Variational :: SafeSigma(double delta) {
	double correction = 0;
	if (Sigma + delta == 0) {
		correction = delta*((double)((int)Rnd->Uniform(0, 2))-0.5)/5.; // Random correction on the left or on the right
	}
	else if (Sigma + delta < 0) {
		return Sigma/2.; 		 // Using delta = -Sigma/2.
	} 
	//cout << "CORRECTION = " <<  correction << endl;
	return Sigma + delta + correction;
}

int Variational :: MoveMean() {
	int C = 1; // Cycle counter
	double mleft, mright; // Left, Right mean
	double hleft, hright; // Left, Right <H>
	double *SampleLeft = new double[Nsample];
	double *SampleRight = new double[Nsample];
	
	// Compute H at (Mean, Sigma)
	GenerateUniform();
	ComputeExpH();
	
	// Compute H at (Left mean, Sigma)
	mleft = Mean - R_Mean;
	GenerateUniform(SampleLeft, mleft, Sigma);
	hleft = ComputeExpH(SampleLeft, mleft, Sigma);
	
	// Compute H at (Right mean, Sigma)
	mright = Mean + R_Mean;
	GenerateUniform(SampleRight, mright, Sigma);
	hright = ComputeExpH(SampleRight, mright, Sigma);
		
	do {
		if (hleft < h) {
			if (hright >= h) { // MOVE LEFT
				//if (Print) PrintMoveMean("LEFT", C, mleft, mright, hleft, hright);
				// SHIFT
				mright = Mean;
				Mean = mleft;
				hright = h;
				h = hleft;
				SampleRight = Sample;
				Sample = SampleLeft; 
				// GENERATE NEW
				mleft = Mean - R_Mean;
				GenerateUniform(SampleLeft, mleft, Sigma);
				hleft = ComputeExpH(SampleLeft, mleft, Sigma);
			}
			else {
				if (hleft < hright) { // MOVE LEFT	
					//if (Print) PrintMoveMean("LEFT", C, mleft, mright, hleft, hright);
					// SHIFT
					mright = Mean;
					Mean = mleft;
					hright = h;
					h = hleft;
					SampleRight = Sample;
					Sample = SampleLeft; 
					// GENERATE NEW
					mleft = Mean - R_Mean;
					GenerateUniform(SampleLeft, mleft, Sigma);
					hleft = ComputeExpH(SampleLeft, mleft, Sigma);
				}
				else { // MOVE RIGHT
					//if (Print) PrintMoveMean("RIGHT", C, mleft, mright, hleft, hright);
					// SHIFT
					mleft = Mean;
					Mean = mright;
					hleft = h;
					h = hright;
					SampleLeft = Sample;
					Sample = SampleRight;
					// GENERATE NEW
					mright = Mean + R_Mean;
					GenerateUniform(SampleRight, mright, Sigma);
					hright = ComputeExpH(SampleRight, mright, Sigma);
				}
			}
		}
		else {
			if (hright < h) { // MOVE RIGHT
				//if (Print) PrintMoveMean("RIGHT", C, mleft, mright, hleft, hright);
				// SHIFT
				mleft = Mean;
				Mean = mright;
				hleft = h;
				h = hright;
				SampleLeft = Sample;
				Sample = SampleRight;
				// GENERATE NEW
				mright = Mean + R_Mean;
				GenerateUniform(SampleRight, mright, Sigma);
				hright = ComputeExpH(SampleRight, mright, Sigma);
			}
			else { // STALEMATE: exit 
				if (Print) PrintMoveMean("RETURN", C, mleft, mright, hleft, hright);
				return -1;
			}
		}
		C++;
	} while (true);
	
	return 0;
}

int Variational :: MoveSigma() {
	int C = 1;
	double sleft, sright;
	double hleft, hright;
	double *SampleLeft = new double[Nsample];
	double *SampleRight = new double[Nsample];
	
	// Sigma, h
	GenerateUniform();
	ComputeExpH();
	
	// LEFT
	sleft = SafeSigma(-R_Sigma);
	GenerateUniform(SampleLeft, Mean, sleft);
	hleft = ComputeExpH(SampleLeft, Mean, sleft);
	
	// RIGHT
	sright = Sigma + R_Sigma;
	GenerateUniform(SampleRight, Mean, sright);
	hright = ComputeExpH(SampleRight, Mean, sright);
	
	do {
		if (hleft < h) {
			if (hright >= h) { // MOVE LEFT
				//if (Print) PrintMoveSigma("LEFT", C, sleft, sright, hleft, hright);
				// SHIFT
				sright = Sigma;
				Sigma = sleft;
				hright = h;
				h = hleft;
				SampleRight = Sample;
				Sample = SampleLeft; 
				// GENERATE NEW
				sleft = SafeSigma(-R_Sigma);
				GenerateUniform(SampleLeft, Mean, sleft);
				hleft = ComputeExpH(SampleLeft, Mean, sleft);
			}
			else {
				if (hleft < hright) { // MOVE LEFT
					//if (Print) PrintMoveSigma("LEFT", C, sleft, sright, hleft, hright);
					// SHIFT
					sright = Sigma;
					Sigma = sleft;
					hright = h;
					h = hleft;
					SampleRight = Sample;
					Sample = SampleLeft; 
					// GENERATE NEW
					sleft = SafeSigma(-R_Sigma);
					GenerateUniform(SampleLeft, Mean, sleft);
					hleft = ComputeExpH(SampleLeft, Mean, sleft);
				}
				else { // MOVE RIGHT
				
					//if (Print) PrintMoveSigma("RIGHT", C, sleft, sright, hleft, hright);
					// SHIFT
					sleft = Sigma;
					Sigma = sright;
					hleft = h;
					h = hright;
					SampleLeft = Sample;
					Sample = SampleRight;
					// GENERATE NEW
					sright = Sigma + R_Sigma;
					GenerateUniform(SampleRight, Mean, sright);
					hright = ComputeExpH(SampleRight, Mean, sright);
				}
			}
		}
		else {
			if (hright < h) { // MOVE RIGHT
				//if (Print) PrintMoveSigma("RIGHT", C, sleft, sright, hleft, hright);
				// SHIFT
				sleft = Sigma;
				Sigma = sright;
				hleft = h;
				h = hright;
				SampleLeft = Sample;
				Sample = SampleRight;
				// GENERATE NEW
				sright = Sigma + R_Sigma;
				GenerateUniform(SampleRight, Mean, sright);
				hright = ComputeExpH(SampleRight, Mean, sright);
			}
			else { // STALEMATE: try again with half step, then twice the step; if it fails again abort.
				if (Print) PrintMoveSigma("RETURN", C, sleft, sright, hleft, hright);
				return -1;
			}
		}
		C++;
	} while (true);
	
	return 0;
	
}

void Variational :: Minimize() {
	cout << "-------- Minimization of the energy --------" << endl;
	if (!OptMinimize) return; // Skip minimization process
	
	int ret1, ret2;
	int Stalemate1 = 0, Stalemate2 = 0;
	double Nsample0 = Nsample;
	ResizeSample(NsampleMinim);
	double Mean0, Sigma0;
	double R_Mean0, R_Sigma0;
	double R_Mean_cycle, R_Sigma_cycle;
	Print = false;
	string option;
	
	Mean0 = Mean;
	Sigma0 = Sigma;
	R_Mean0 = R_Mean;
	R_Sigma0 = R_Sigma;
	
	
	
	do {
		
		R_Mean_cycle = R_Mean;   	// Starting point of each total cycle
		R_Sigma_cycle = R_Sigma;
		for (unsigned int i = 0; i < TotalCycles; i++) {
			Stalemate1 = 0;
			Stalemate2 = 0;
			R_Mean = R_Mean_cycle;
			R_Sigma = R_Sigma_cycle;
			for (unsigned int k = 0; k < MoveMeanCycles; k++) {
				//cout << "--------------------------------------------------------------------------" << endl << endl;
				ret1 = MoveMean(); ret1 *= 1; 
				Stalemate1 += ret1;
				R_Mean /= MoveMeanScale;
			} 
			for (unsigned int k = 0; k < MoveSigmaCycles; k++) {
				ret2 = MoveSigma(); ret2 *= 1; 
				Stalemate2 += ret2;
				R_Sigma /= MoveSigmaScale;
				//cout << "--------------------------------------------------------------------------" << endl << endl;
			}
			if (Stalemate1 < -(double)MoveMeanCycles/2.) {
				Mean += Rnd->Uniform(-R_Mean/5., R_Mean/5.);
				GenerateUniform();
				ComputeExpH();
			}
			if (Stalemate2 < -(double)MoveSigmaCycles/2.) {
				double delta = Rnd->Uniform(-R_Sigma/5., R_Sigma/5.);
				SafeSigma(delta);
				GenerateUniform();
				ComputeExpH();
			}
			
			R_Mean_cycle /= CycleScale;
			R_Sigma_cycle /= CycleScale;
			
			cout << "-------- [Cycle " << i+1 << "] Parameters --------" << endl;
			cout << "  \u03BC = " << left << setw(8) << Mean << endl;
			cout << "  \u03C3 = " << left << setw(8) << Sigma << endl;
			cout << "  <H> = " << h << endl;
			cout << endl;
		}
		
		
		cout << "Accept these parameters? (y|n) ";
		cin >> option;
		cout << endl;
		for_each(option.begin(), option.end(), [](char & c){ c = tolower(c);} );
		
		if (option != "y") { // Return to the previous setup
			Mean = Mean0;
			Sigma = Sigma0;
			R_Mean = R_Mean0;
			R_Sigma = R_Sigma0;
		}
		else { // If they are accepted, scale the move step to help the convergence
			R_Mean0 /= CycleScale;
			R_Sigma0 /= CycleScale;
			R_Mean = R_Mean0;
			R_Sigma = R_Sigma0;
		}
		
		cout << "Start another cycle? (y|n) ";
		cin >> option;
		cout << endl;
		for_each(option.begin(), option.end(), [](char & c){ c = tolower(c);} );
		if (option != "y") break;
		else {
			
		}
	} while(true);
	cout << endl;
	//R_Mean = R_Mean0;
	//R_Sigma = R_Sigma0;
	ResizeSample(Nsample0);
	GenerateUniform();
	ComputeExpH();
	return;
}

bool Variational :: FileExists(const string& filename) {
	struct stat buffer;   
	return (stat(filename.c_str(), &buffer) == 0); 
}

void Variational :: Input(const string filename) {
	if (!FileExists(filename)) {
		cerr << "Error [Variational :: Input(const string)]: cannot open the input file" << endl;
		exit(1);
	}
	
	ifstream In(filename);
	string Option;
	//double r;
	int n, l;
	In.ignore(250, '\n');
	In >> n;								In.ignore(250, '\n'); // Sample size (number of points)
	ResizeSample(n);
	In >> n;								In.ignore(250, '\n'); // Number of blocks
	In >> l;								In.ignore(250, '\n'); // Block Length
	ResizeStat(n, l);
	In >> Mean; 							In.ignore(250, '\n');
	In >> Sigma; 							In.ignore(250, '\n');
	if (Sigma <= 0) {						In.ignore(250, '\n');
		cerr << "Error [Variational :: Input(const string)]: invalid Sigma" << endl;
		exit(1);
	}
	In >> Start; 							In.ignore(250, '\n'); // Starting point
	In >> R_Psi; 							In.ignore(250, '\n');
	if (R_Psi <= 0) {						In.ignore(250, '\n'); // Metropolis range of motion
		cerr << "Error [Variational :: Input(const string)]: invalid R_Psi" << endl;
		exit(1);
	}
	
	//------ MINIMIZATION PARAMETERS ------//
	In.ignore(250, '\n');
	In >> Option;							In.ignore(250, '\n');
	for_each(Option.begin(), Option.end(), [](char & c){ c = tolower(c);} );
	if (Option == "true") OptMinimize = true;
	else OptMinimize = false;
	
	In >> NsampleMinim; 					In.ignore(250, '\n'); // Sample size (Minimize() function)
	if (NsampleMinim <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid NsampleMinim" << endl;
		exit(1);
	}
	In >> TotalCycles;						In.ignore(250, '\n'); 
	if (TotalCycles <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid TotalCycles" << endl;
		exit(1);
	}
	In >> MoveMeanCycles; 					In.ignore(250, '\n');
	if (MoveMeanCycles <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid MoveMeanCycles" << endl;
		exit(1);
	}
	In >> MoveSigmaCycles; 					In.ignore(250, '\n');
	if (MoveSigmaCycles <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid MoveSigmaCycles" << endl;
		exit(1);
	}
	In >> R_Mean;   						In.ignore(250, '\n');
	if (R_Mean <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid R_Mean" << endl;
		exit(1);
	}
	In >> R_Sigma; 							In.ignore(250, '\n');  
	if (R_Sigma <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid R_Sigma" << endl;
		exit(1);
	}
	In >> CycleScale; 						In.ignore(250, '\n');
	if (CycleScale <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid CycleScale" << endl;
		exit(1);
	}
	In >> MoveMeanScale; 					In.ignore(250, '\n');
	if (MoveMeanScale <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid MoveMeanScale" << endl;
		exit(1);
	}
	In >> MoveSigmaScale; 					In.ignore(250, '\n');		
	if (MoveSigmaScale <= 0) {
		cerr << "Error [Variational :: Input(const string)]: invalid MoveSigmaScale" << endl;
		exit(1);
	}
	In.close();
}

void Variational :: PrintMoveMean(string s, int C, double mleft, double mright, double hleft, double hright) {
	int sp = 8;
	cout << fixed << setprecision(sp-2);
	cout << "[Cycle = " << C << "]   ";
	cout << s <<"   MEAN = [" <<setw(sp)<< mleft << " | " <<setw(sp)<< Mean << " | " <<setw(sp)<< mright << "] ----> [" <<setw(sp)<< hleft << " | " <<setw(sp)<< h << " | " <<setw(sp)<< hright << "]" << endl;
	return ;
}

void Variational :: PrintMoveSigma(string s, int C, double sleft, double sright, double hleft, double hright) {
	int sp = 8;
	cout << fixed << setprecision(sp-2);
	cout << "[Cycle = " << C << "]   ";
	cout << s << "   SIGMA = [" <<setw(sp)<< sleft << " | " <<setw(sp)<< Sigma << " | "<<setw(sp)<< sright << "] ----> [" <<setw(sp)<< hleft << " | " <<setw(sp)<< h << " | " <<setw(sp)<< hright << "]" << endl;
	return ;
}


void Variational :: Simulate() {
	cout << "-------- Simulation --------" << endl;
	for (unsigned int block = 0; block < Nblocks; block++) {
		cout << "[Block " << block+1 << "] " << endl;
		for (unsigned int step = 0; step < Lblock; step++) {
			GenerateUniform();
			ComputeExpH();
			H.Append(h);
		}
	}
	cout << endl;
	vector<double> prog_average, prog_sigma;
	H.BlockProgAverageSigma(prog_average, prog_sigma);
	OutputProgH(prog_average, prog_sigma);
	return;
}

void Variational :: Calibrate() {
	cout << "-------- Calibration of the acceptance ratio --------" << endl;
	char input[10];
	ResizeSample(Nsample);
	GenerateUniform();
	cout.precision(4);
	double input_radius;
	do {
		cout << "    Actual range of motion = " << fixed << R_Psi << ", Acceptance = " << RatioAcc*100. << "%"
			 << endl << "    Insert the new radius or quit calibration (q): ";
		cin >> input;
		cout << endl;
		try { input_radius = stod(input); }
		catch(std::invalid_argument) { break; }
		
		if (input_radius > 0) R_Psi = input_radius;
		GenerateUniform();
	} while(true);
	return;
}

void Variational :: AppendH() {
	H.Append(h);
	return;
}

void Variational :: OutputProgH(vector<double> prog_average, vector<double> prog_sigma) const {
	cout << "-------- Saving the progressive value of <H> on <OutputData/H_prog.out>" << endl;
	ofstream Out("OutputData/H_prog.out");
	Out << "Nblocks: " << Nblocks << endl;
	Out << "BlockLength: " << Lblock << endl;
	Out << "Mean: " << Mean << endl;
	Out << "Sigma: " << Sigma << endl;
	for (unsigned int block = 0; block < Nblocks; block++) {
		Out << prog_average[block] << "    " << prog_sigma[block] << endl;
	}
	Out.close();
	return;
}


void Variational :: OutputSample() const {
	cout << "-------- Saving the sample on <OutputData/sample.out>" << endl;
	ofstream Out("OutputData/sample.out");
	for (unsigned int i = 0; i < Nsample; i++) {
		Out << Sample[i] << endl;
	}
	Out.close();
}

void Variational :: PrintSample() const {
	for (unsigned int i = 0; i < Nsample; i++) {
		cout << "[" << i << "] = " << Sample[i] << endl;
	}
}

void Variational :: PrintSample(double* SampleNew) const {
	for (unsigned int i = 0; i < Nsample; i++) {
		cout << "[" << i << "] = " << SampleNew[i] << endl;
	}
}

