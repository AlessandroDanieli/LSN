#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>

#include "random.h"
#include "stat.h"
#include "metropolisorbitals.h"
#include "auxiliary.h"
using namespace std;

int main(int argc, char** argv) {
	Random Rnd;		Rnd.RandGenSetup();
	const double a0 = 0.0529; // Bohr Radius in nanometers
	
	unsigned int Nsampling = 2; // Uniform, Gaussian.
	string FilenameMethod[Nsampling] = {"_unif", "_gauss"};
	
	int **QN, ret;
	unsigned int Norbitals, Nsample, Nblocks, CalibSteps, Sample_first, Sample_last;
	Vector* Start;
	double *RangeOfMotion, r_mean;
	bool *CalibOption, ExistCalibTrue = false;
	
	// Read input parameters
	ret = Input("input.dat", QN, Norbitals, Nsample,Nblocks,Start,RangeOfMotion,CalibOption,CalibSteps,Nsampling,Sample_first,Sample_last);
	if (ret == -1) {
		cerr << "Error: cannot open the file <input.dat>" << endl;
		exit(-1);
	}
	cout << "Sample size: " << Nsample << endl;
	cout << "Number of blocks: " << Nblocks << endl;
	cout << "Number of orbitals: " << Norbitals << endl << endl;
	
	double TimeCalib[Nsampling*Norbitals], TimeComp[Nsampling*Norbitals];
	string folder = "OutputData/", filename; // Output folder, filename
	vector<double> prog_average, prog_sigma;
	Stat Radius(Nsample, Nblocks);
	MetropolisOrbital M(&Rnd, 3, Nsample); // Metropolis Algorithm	
	
	// Calibration cycles
	for (unsigned int n = 0; n < Norbitals; n++) {
		M.SetQN(QN[n][0], QN[n][1], QN[n][2], QN[n][3]);
		M.SetStart(Start[n]);
		if (CalibOption[n]) {
			if (!ExistCalibTrue) cout << "****** Calibration ******" << endl;
			ExistCalibTrue = true;
			cout << "[" << n+1 << "] WaveFunction \u03A8[" << QN[n][0] << ", "<<QN[n][1]<<", "<<QN[n][2]<<", "<<QN[n][3]<<"]" << endl; 
			
			for (unsigned int sampling = 0; sampling < Nsampling; sampling++) {
				TimeCalib[Nsampling*n+sampling] = M.Calibrate(sampling, RangeOfMotion[Nsampling*n + sampling], CalibSteps);
			}			
		}
	}
	
	cout << endl << "****** Computation ******" << endl;
	for (unsigned int n = 0; n < Norbitals; n++) {
		M.SetQN(QN[n][0], QN[n][1], QN[n][2], QN[n][3]);
		M.SetStart(Start[n]);	
		filename = folder+to_string(QN[n][0])+", "+to_string(QN[n][1])+", "+to_string(QN[n][2])+", "+to_string(QN[n][3]);
		cout << "[" << n+1 << "] WaveFunction \u03A8[" << QN[n][0] << ", "<<QN[n][1]<<", "<<QN[n][2]<<", "<<QN[n][3]<<"]" << endl; 
		
		for (unsigned int sampling = 0; sampling < Nsampling; sampling++) {
			TimeComp[Nsampling*n+sampling] = M.Generate(sampling, RangeOfMotion[Nsampling*n + sampling]);
			r_mean = M.MeanRadius();
			cout << setprecision(4) << fixed;
			cout << "    Mean Radius: " << r_mean << " * a0 = " << r_mean * a0 << " nm" << endl;
			cout << "    Acceptance ratio: " << M.GetRatioAcc()*100 << " %" << endl;
			cout << "    Analysing the radial distribution..." << endl;
			
			Radius.CopyData(M.GetRadius());
			Radius.BlockProgAverageSigma(prog_average, prog_sigma);
			cout << "    Outputting the data..." << endl << endl;
			filename = folder + to_string(QN[n][0])+"."+to_string(QN[n][1])+"."+to_string(QN[n][2])+"."+to_string(QN[n][3]);
			filename = (filename+FilenameMethod[sampling]).c_str();
			OutputProgRadius(prog_average, prog_sigma, filename+".out", Nsample);
			M.OutputSample(filename+".xyz", Sample_first, Sample_last);
		}
	}
	
	// Update the input file with the new "Range of motion" radii. Even if these radii are left untouched, it might be useful 
	// to align again the data in a table-like structure
	OutputFinalSetup("input.dat", QN, Norbitals, Nsample, Nblocks, Start, RangeOfMotion, CalibOption, CalibSteps, Sample_first, Sample_last);

	Rnd.SaveSeed();
	return 0;
}
