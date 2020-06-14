#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "stat.h"
#include "auxiliary.h"
using namespace std;
 
int main (int argc, char *argv[]) {
	Random rnd;
	rnd.RandGenSetup();
	
	//## Parameters ##//
	unsigned int M = 10000000; // Number of measures
	unsigned int N = 250; // Blocks
	unsigned int L = M/N;
	vector <double> prog_average, prog_error, chi_array;
	Stat s(&rnd, M, N);
	double t1 = s.FillUniform();
	
	//## Data Analysis - AVERAGE ##//
	string filename = "average.out";
	double t2 = s.BlockProgAverage(prog_average, prog_error);
	double ot1 = ProgOutput(filename, prog_average, prog_error, L);
	
	//## Data Analysis - RMS ##//
	filename = "error.out";
	s.SetExpMean(0.5);
	double t3 = s.BlockProgSigmaExp(prog_average, prog_error);
	double ot2 = ProgOutput(filename, prog_average, prog_error, L);
	
	//## Data Analysis - CHI test ##//
	M = 100; // Chi sample: 100 blocks
	L = 10000; // Length of each block
	double chi2;
	//unsigned int bins = 20; // Chi2 test bins
	double t4 = Chi(s, chi_array, chi2, M, L);
	double ot3 = ChiOutput("chi.out", chi_array, chi2, M, L);
	
	cout << "Elapsed Times: " << endl;
	cout << "> FillUniform = " << t1 << " \u03BCs" << endl;
	cout << "> BlockProgAverage = " << t2 << " \u03BCs" << endl;
	cout << "> ProgProgSigmaExp = " << t3 << " \u03BCs" << endl;
	cout << "> Chi = " << t4 << " \u03BCs" << endl;
	
	cout << "> ProgOutput = " << ot1 << " \u03BCs" << endl;
	cout << "> ProgOutput = " << ot2 << " \u03BCs" << endl;
	cout << "> ChiOutput = " << ot3 << " \u03BCs" << endl;
	rnd.SaveSeed();
	return 0;
}

