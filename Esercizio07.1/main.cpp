#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"
#include "stat.h"
#include "canonical.h"
using namespace std;

void Input(const string inputFile, Canonical& C) {
	ifstream Input;
	// INPUT
	unsigned int n1, n2;
	double d, rcut, rho, vtail, ptail;
	Input.open(inputFile);
	Input.ignore(120, '\n');
	Input >> n1; 		Input.ignore(120, '\n'); 	C.SetNp(n1);
	Input >> n1; 		Input.ignore(120, '\n');
	Input >> n2;		Input.ignore(120, '\n');	C.SetBlocksSteps(n1, n2);
	Input.ignore(120, '\n');
	Input >> d; 		Input.ignore(120, '\n'); 	C.SetT(d);
	Input >> rho;		Input.ignore(120, '\n'); 	C.SetRho(rho); 	C.SetBox(); // Sets volume, L_box
	Input >> rcut;		Input.ignore(120, '\n'); 	C.SetRcut(rcut);
	Input >> d;			Input.ignore(120, '\n'); 	C.SetRangeOfMotion(d);
	Input.close();
	
	cout << "Reading input parameters from <" << inputFile << ">:" << endl;
	cout << "> Number of particles = " << C.GetNp() << endl;
	cout << "> Number of blocks = " << C.GetNb() << endl;
	cout << "> Number of steps for each block = " << C.GetNs() << endl;
	//cout << "> Number of bins (g(r) sampling)= " << C.GetNbins() << endl;
	cout << "> Temperature = " << C.GetT() << endl;
	cout << "> Density of particles = " << C.GetRho() << endl;
	cout << "> Volume of the simulation box = " << C.GetVol() << endl;
	cout << "> Edge of the simulation box L = " << C.GetL() << endl;
	cout << "> Cutoff radius of the interatomic potential = " << C.GetRcut() << endl << endl;

	//Tail corrections for potential energy and pressure
	vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
	cout << "> Tail correction for the potential energy = " << vtail << endl;
	cout << "> Tail correction for the virial           = " << ptail << endl; 
	cout << "> Using Metropolis algorithm (Range of motion = " << C.GetRangeOfMotion() << ")" << endl << endl;
	return;
}


int main() {
	Random Rnd; Rnd.RandGenSetup();
	Canonical C(&Rnd);
	cout << "Classic Lennard-Jones fluid (using Lennard-Jones units)" << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
	cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl;
	
	Input("input.dat", C);
	C.LoadConfig("config.0");
	C.PrintParameters();
	
	double RatioAcc, InputRange;
	string Input;
	do {
		RatioAcc = C.Calibration(300);
		cout << "    Actual range of motion = " << C.GetRangeOfMotion() << ", Acceptance = " << RatioAcc*100. << "%"
			 << endl << "    Insert the new radius or quit calibration (q): ";
		cin >> Input;
		cout << endl;
		try { InputRange = stod(Input); }
		catch(std::invalid_argument) { break; }
		if (InputRange > 0) C.SetRangeOfMotion(InputRange);
	} while(true);
	
	C.AutoEquilibration(400, 0.012);
	C.Simulate("OutputData/inst_UP.out");

	Rnd.SaveSeed();
	return 0;
}
