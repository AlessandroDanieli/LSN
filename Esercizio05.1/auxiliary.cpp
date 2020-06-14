#include "auxiliary.h"
#include <iomanip>
using namespace std;

bool FileExists(const string& filename) {
	struct stat buffer;   
	return (stat(filename.c_str(), &buffer) == 0); 
}


int Input(string filename, int**& QN, unsigned int& Norbitals, unsigned int& Nsample, unsigned int& Nblocks, Vector*& Start, double*& RangeOfMotion, bool*& CalibOption, unsigned int& CalibSteps, unsigned int Nsampling, unsigned int& Sample_first, unsigned int& Sample_last){
	if (!FileExists(filename)) {
		return -1;
	}
	else {
		ifstream Input(filename);
		if (Input.is_open()) {
			string Option, buff;
			double x, y, z;
			Input >> Nsample;						Input.ignore(250, '\n'); // Sample size (number of points)
			Input >> Nblocks;						Input.ignore(250, '\n'); // Number of blocks
			Input >> Sample_first >> Sample_last; 	Input.ignore(250, '\n'); // Range of the sample to be saved on file
			Input >> CalibSteps;					Input.ignore(250, '\n'); // Calibration steps
			
			// Orbitals
			Input.ignore(250, '\n');
			Input.ignore(250, '\n');
			Input >> Norbitals;						Input.ignore(250, '\n'); // Number of orbitals
			Input.ignore(250, '\n');
			Input.ignore(250, '\n');
			
			QN = new int*[Norbitals];
			RangeOfMotion = new double[Nsampling*Norbitals];
			CalibOption = new bool[Norbitals];
			Start = new Vector[Norbitals];
			
			for (unsigned int i = 0; i < Norbitals; i++) {
				Start[i].SetDim(3);
				QN[i] = new int[4];
				for (unsigned int k = 0; k < 4; k++) {
					Input >> QN[i][k] >> buff;
				}
				for (unsigned int sampling = 0; sampling < Nsampling; sampling++) {
					Input >> RangeOfMotion[Nsampling*i+sampling] >> buff;
				}
				// Calibration option [true|false]
				Input >> Option >> buff;
				for_each(Option.begin(), Option.end(), [](char & c){ c = tolower(c);} );
				if (Option == "true") CalibOption[i] = true;
				else CalibOption[i] = false;
				
				Input >> x >> buff >> y >> buff >> z;	Input.ignore(250, '\n'); // Starting point coordinates
				Start[i].Set({x, y, z});
			}
			Input.close();
		}
		else {
			cerr << "Error: cannot open the input file." << endl;
			exit(-1);
		}
		
		return 0;
	}
}

/*
int InputQN(const string filename, int**& QN, unsigned int& Norbitals) {
	if (!FileExists(filename)) {
		return -1;
	}
	else {
		ifstream Input(filename);
		Input >> Norbitals;			Input.ignore(250, '\n'); // Number of orbitals
		QN = new int*[Norbitals];
		for (unsigned int i = 0; i < Norbitals; i++) {
			QN[i] = new int[4];
			for (unsigned int k = 0; k < 4; k++) {
				Input >> QN[i][k];
			}
		}
		return 0;
	}
}
*/

// Overwrites the input file subsituting the "Range of motion" values for each orbital achieved by the calibration
void OutputFinalSetup(string filename, int** QN, unsigned int Norbitals, unsigned int Nsample, unsigned int Nblocks, Vector* Start, double* R, bool* CalibOption, unsigned int CalibSteps, unsigned int Sample_first, unsigned int Sample_last) {
	ofstream Out(filename);
	int spacing = 27;
	if (Out.is_open()) {
		Out << left;
		Out << setw(spacing) << Nsample << "# Sample size" << endl;
		Out << setw(spacing) << Nblocks << "# Number of blocks (blocking method)" << endl;
		Out << setw(12) << Sample_first << " " << setw(12) << Sample_last << "  ";
		Out << "# Range of the sample to be saved on file" << endl;
		Out << setw(spacing) << CalibSteps << "# Calibration steps (if the corresponding boolean value is true, see below)" << endl;
		
		Out << "_______________________________________________________________________________________________________________________________" << endl;
		Out << "######  Orbitals  ######" << endl;
		Out << setw(spacing) << Norbitals << "# Number of orbitals to be computed" << endl;
		Out << "_______________________________________________________________________________________________________________________________"<<endl;
		Out << 
"Z    | n    | l    | m    | Uniform radius    | Gaussian radius    | Calibration [true|false] | X start   | Y start   | Z start" << endl;
		
		/*Out << "###     " */
		for (unsigned int n = 0; n < Norbitals; n++) {
			Out << left;
			Out <<setw(5) << QN[n][0] <<"| "<<setw(5) << QN[n][1] <<"| "<<setw(5) << QN[n][2] <<"| "<<setw(5) << QN[n][3] <<"| ";
			// THIS SHOULD BE GENERALIZED (IF Nsampling != 2)!
			Out <<setw(18) << R[2*n] <<"| "<<setw(19) << R[2*n+1] <<"| ";
			if (CalibOption[n]) Out << setw(25) << "true";
			else Out << setw(25) << "false";
			
			Out << "| " << setw(10) << Start[n][0] << "| " << setw(10) << Start[n][1] << "| " << setw(10) << Start[n][2] << endl;
		}
		Out.close();
	}
	else {
		cerr << "Error: cannot open the input file." << endl;
		exit(-1);
	}
}


void OutputProgRadius(vector<double> prog_average, vector<double> prog_sigma, string filename, unsigned int Nsample) {
	ofstream Out(filename);
	if (Out.is_open()) {
		Out << "#Samples " << Nsample << endl;
		Out << "#Blocks " << prog_average.size() << endl;
		for (unsigned int i = 0; i < prog_average.size(); i++) Out << prog_average[i] << " " << prog_sigma[i] << endl;
	}
	else {
		cerr << "Error: cannot open output file." << endl;
		exit(-1);
	}
	Out.close();
	return;
}


void PrintQN(int**& QN, unsigned int& Norbitals) {
	for (unsigned int n = 0; n < Norbitals; n++) {
		cout << QN[n][0] << " " << QN[n][1] << " " << QN[n][2] << " " << QN[n][3] << endl;
	}
	cout << endl;
	return;
}


void PrintRange(double* RangeOfMotion, unsigned int Norbitals) {
	for (unsigned int n = 0; n < Norbitals; n++) {
		cout << RangeOfMotion[2*n] << " " << RangeOfMotion[2*n+1] << endl;
	}
	cout << endl;
	return;
}
