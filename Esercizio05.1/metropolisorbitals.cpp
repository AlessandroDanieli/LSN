#include "metropolisorbitals.h"
using namespace std;

// *** Advice: when plotting it might be useful to introduce a scaling factor for rho (usually /integer) in R_nn2 or other functions *** /

MetropolisOrbital :: MetropolisOrbital(): Orbital(), Rnd(nullptr), Sample(nullptr), d(1), Nsample(0) {}

MetropolisOrbital :: MetropolisOrbital(Random* _Rnd, unsigned int _d, unsigned int _Ns): Orbital(), Rnd(_Rnd), d(_d), Nsample(_Ns) {
	Sample = new Vector[Nsample];
}

// Default T ~ Uniform distribution (symmetric)
MetropolisOrbital :: MetropolisOrbital(Random* _Rnd, unsigned int _Nsample, Vector start): Orbital(), Rnd(_Rnd), Nsample(_Nsample) {
	Sample = new Vector[Nsample];
	Sample[0] = start;
	d = start.GetDim(); // Dimension 
}

MetropolisOrbital :: MetropolisOrbital(unsigned int _Z, unsigned int _n, unsigned int _l, int _m): Orbital(_Z, _n, _l, _m)  {
	MetropolisOrbital();
}

MetropolisOrbital :: MetropolisOrbital(Random* _Rnd, unsigned int _Nsample, Vector start, unsigned int _Z, unsigned int _n, 
									   unsigned int _l, int _m) {					   
	Orbital(_Z, _n, _l, _m);
	MetropolisOrbital(_Rnd, _Nsample, start);
}


MetropolisOrbital :: ~MetropolisOrbital() {}

void MetropolisOrbital :: SetStart(Vector& Start)  {
	d = Start.GetDim();
	Sample[0] = Start;
}

vector<Vector> MetropolisOrbital :: GetSample() const {
	vector<Vector> _Sample;
	for (unsigned int i = 0; i < Nsample; i++) _Sample.push_back(Sample[i]);
	return _Sample;
}

double MetropolisOrbital :: Generate(unsigned int sampling, double& RangeOfMotion) {
	if (sampling == 0) return GenerateUniform(RangeOfMotion);
	else if (sampling == 1) return GenerateGaussian(RangeOfMotion);
	else {
		cerr << "Error [MetropolisOrbital :: Generate(unsigned int, double&)]: Invalid sampling method." << endl;
		exit(-1);
	}
}

double MetropolisOrbital :: Calibrate(unsigned int sampling, double& RangeOfMotion, unsigned int CalibSteps) {
	if (sampling == 0) return CalibrateUniform(RangeOfMotion, CalibSteps);
	else if (sampling == 1) return CalibrateGaussian(RangeOfMotion, CalibSteps);
	else {
		cerr << "Error [MetropolisOrbital :: Calibrate(unsigned int, double&, unsigned int)]: Invalid sampling method." << endl;
		exit(-1);
	}
}


// ## Uniform Sampling (CUBIC) ## //
double MetropolisOrbital :: GenerateUniform(double& RangeOfMotion) {
	auto t1 = chrono::high_resolution_clock::now();
	unsigned int Naccepted = 0;
	double q;
	Vector X_try(d);
	
	for (unsigned int i = 1; i < Nsample; i++) {
		for (unsigned int k = 0; k < d; k++) {
			X_try[k] = Sample[i-1][k] + Rnd->Uniform(-RangeOfMotion, RangeOfMotion);
		}				
		
		q = pv_nn(X_try)/pv_nn(Sample[i-1]);
		
		if (Rnd->Uniform() < q) {
			Sample[i].Copy(X_try); 
			Naccepted++;	
		}
		else {
			Sample[i].Copy(Sample[i-1]);
		}	
	}
	RatioAcc = Naccepted/(double)(Nsample-1);
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}



double MetropolisOrbital :: GenerateGaussian(double& RangeOfMotion) {
	auto t1 = chrono::high_resolution_clock::now();
	unsigned int Naccepted = 0;
	double q; // Ratio p(x_try)/p(x_prev)
	Vector X_try(d);
	
	for (unsigned int i = 1; i < Nsample; i++) {
		for (unsigned int k = 0; k < d; k++) {
			X_try[k] = Sample[i-1][k] + Rnd->Gauss(0, RangeOfMotion); // Sigma is a variable
		}
		q = pv_nn(X_try)/pv_nn(Sample[i-1]);
		
		if (Rnd->Uniform() < q) {
			Sample[i].Copy(X_try); 
			Naccepted++;
		}
		else {
			Sample[i].Copy(Sample[i-1]);
		}
	}
	RatioAcc = Naccepted/(double)(Nsample-1);
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}

double MetropolisOrbital :: CalibrateUniform(double& RangeOfMotion, unsigned int CalibSteps) { // Reaches a 50% acceptance ratio with <delta> as tolerance
	auto t1 = chrono::high_resolution_clock::now();
	
	// Temporary substitution to limit the number of iterations in the calibration process (using <Ncalib> steps)
	char input[10];
	unsigned int MemNsample = Nsample;
	
	cout << "    [Uniform sampling]" << endl;
	// Temporarily set Nsample = Ncalib in order to move via GenerateUniform <Nsample> times
	if (CalibSteps < 10000) Nsample = 10000; // Minimum number of iterations
	if (CalibSteps > Nsample) Nsample = MemNsample; // Avoid exceeding Sample's size
	else Nsample = CalibSteps;

	GenerateUniform(RangeOfMotion);
	cout.precision(4);
	double input_radius;
	do {
		cout << "    Actual range of motion = " << RangeOfMotion << ", Acceptance = " << RatioAcc*100. << "%"
			 << endl << "    Insert the new radius or quit calibration (q): ";
		cin >> input;
		cout << endl;
		try { input_radius = stod(input); }
		catch(std::invalid_argument) { break; }
		
		if (input_radius > 0) RangeOfMotion = input_radius;
		GenerateUniform(RangeOfMotion);
	} while(true);
	// Restore the original value of <Nsample>
	Nsample = MemNsample;
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}

double MetropolisOrbital :: CalibrateGaussian(double& RangeOfMotion, unsigned int CalibSteps) { // Reaches a 50% acceptance ratio with <delta> as tolerance
	auto t1 = chrono::high_resolution_clock::now();
	// Temporary substitution to limit the number of iterations in the calibration process (using <Ncalib> steps)
	char input[10];
	unsigned int MemNsample = Nsample;
	
	cout << "    [Gaussian sampling]" << endl;
	// Temporarily set Nsample = Ncalib in order to move via GenerateGaussian <Nsample> times
	if (CalibSteps < 10000) Nsample = 10000; // Minimum number of iterations
	if (CalibSteps > Nsample) Nsample = MemNsample; // Avoid exceeding Sample's size
	else Nsample = CalibSteps;
	
	GenerateGaussian(RangeOfMotion);
	cout.precision(4);
	double input_radius;
	do {
		cout << "    Actual range of motion = " << RangeOfMotion << ", Acceptance = " << RatioAcc*100. << "%"
			 << endl << "    Insert the new radius or quit calibration (q): ";
		cin >> input;
		cout << endl;
		try {input_radius = stod(input);}
		catch(std::invalid_argument) {break;}
		
		if (input_radius > 0) RangeOfMotion = input_radius;
		
		GenerateGaussian(RangeOfMotion);
	} while(true);
	// Restore the original value of <Nsample>
	Nsample = MemNsample;
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
}



vector<double> MetropolisOrbital :: GetRadius() const {
	vector<double> R;
	for (unsigned int i = 0; i < Nsample; i++) {
		R.push_back(Sample[i].Norm());
	}
	return R;
}

double MetropolisOrbital :: MeanRadius() const {
	double r = 0;
	for (unsigned int i = 0; i < Nsample; i++) r += Sample[i].Norm();
	return r / (double)Nsample;
}

void MetropolisOrbital :: Print() const {
	for (unsigned int i = 0; i < Nsample; i++) {
		cout << "[" << i+1 << "] -> ";
		Sample[i].Print();
	}
	cout << endl;
}

void MetropolisOrbital :: OutputSample(string filename, unsigned int first, unsigned int last, unsigned int step) const {
	if (first > last || first >= Nsample || last > Nsample || step >= abs((int)(last-first))) {
		cerr << "Error [ Metropolis :: OutputSample(string, unsigned int, unsigned int, unsigned int)]: invalid function parameters"<< endl;
		exit(-1);
	}
	ofstream Out(filename);
	unsigned int res = (last-first)%step;
	if (res != 0) last -= res;
	if (Out.is_open()) {
		for (unsigned int i = first; i < last; i+=step) {
			for (unsigned int k = 0; k < d-1; k++) Out << Sample[i][k] << " ";
			Out << Sample[i][d-1] << endl;
		}
	}
	else {
		cerr << "Error [ Metropolis :: OutputSample(string, unsigned int, unsigned int, unsigned int)]: cannot open the output file"<< endl;
		exit(-1);
	}
	Out.close();
	return;
}
