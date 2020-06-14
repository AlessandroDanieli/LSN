#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include "moldyn.h"
using namespace std;

int Input(const string filename, MolDyn& M, bool& old_config, unsigned int& NeqCycles, unsigned int& NeqSteps) {
	ifstream Input;
	Input.open(filename);
	unsigned int Np, Ns, Nframe, Nmeasure, Nprint, Nblocks;
	double T, rho, rcut, delta;
	string _restart, _old_config;
	if (Input.is_open()) {
		Input.ignore(250, '\n');
		Input >> Np;			Input.ignore(250, '\n'); // Input number of particles
		Input >> Ns;  			Input.ignore(250, '\n'); // Number of steps
		Input >> Nframe; 		Input.ignore(250, '\n'); // Output frequency of frames
		Input >> Nmeasure; 		Input.ignore(250, '\n'); // Output frequency of measures
		Input >> Nblocks;		Input.ignore(250, '\n'); // Number of blocks (Data blocking for the output measures)
		Input >> Nprint; 		Input.ignore(250, '\n'); // Current step info frequency
		Input.ignore(250, '\n');
		Input >> NeqCycles;		Input.ignore(250, '\n'); // Number of equilibration cycles
		Input >> NeqSteps;		Input.ignore(250, '\n'); // Number of steps for each equilibration cycle
		Input.ignore(250, '\n');
		Input >> T;  			Input.ignore(250, '\n'); // Input temperature
		Input >> rho; 			Input.ignore(250, '\n'); // Input volume
		Input >> rcut;  		Input.ignore(250, '\n'); // Cutoff radius
		Input >> delta;  		Input.ignore(250, '\n'); // Time interval
		Input.ignore(250, '\n');
		
		if ((int)((double)Ns/(double)(Nmeasure))%(int)Nblocks != 0) {
			cerr << "Warning: the number of blocks is not a divisor of the number of measures. The simulation is aborted." << endl;
			return -1;
		}
		
		Input >> _old_config; 	Input.ignore(120, '\n');
		std::for_each(_old_config.begin(), _old_config.end(), [](char & c){ c = tolower(c);} );
		if (_old_config == "true") old_config = true;
		else old_config = false;
		/*
		Input >> _restart; 		Input.ignore(120, '\n');
		std::for_each(_restart.begin(), _restart.end(), [](char & c){ c = tolower(c);} );
		if (_restart == "true") restart = true;
		else restart = false;
		*/
		
		M.SetParameters(Np, T, rho, rcut, delta, Ns, Nframe, Nmeasure, Nprint, Nblocks);
		Input.close();
 	}
 	else {
		cerr << "Error [Main :: Input]: cannot open the file <" << filename << ">" << endl;
		exit(-1);
	}
	return 0;
}

int main(int argc, char** argv) {
	cout << "> Classic Lennard-Jones fluid - NVE ensemble - Verlet method" << endl;
	cout << "> Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6] (LJ units)" << endl << endl;
	
	bool old_config;
	string option; 
	int r, r1, r2;
	double t1, t2, t3; // Elapsed times
	unsigned int NeqCycles, NeqSteps;
	MolDyn M; 
	
	
	// Program cycles. Allows to perform multiple simulations, by changin only the input file "input.dat"
	for (unsigned int cycle = 0;  ; cycle++) {
		// Read the file "input.dat" and set all the fundamental parameters
		r = Input("input.dat", M, old_config, NeqCycles, NeqSteps);
		
		// Clear the old frames, if the simulation is going to produce new ones
		if (M.GetOutputF()) r = system("rm -rf frames/*.xyz");
		
		if (r == 0) { // Correct number of blocks (divisor of the number of measures). The calculation can be performed
			if (old_config == true) {
				r1 = M.SetConfig("old.final");
				r2 = M.SetOldConfig("old.0");
				
				if (r1 != 0) {
					cerr << "Warning: configuration file old.final not available" <<endl;
					r1 = M.SetConfig("config.0");
					if (r1 != 0) {
						cerr << "Error: configuration file config.0 not available. Aborting the simulation" <<endl;
						return -1;
					}
					cerr << "Warning: cannot read the file old.final. Using instead config.0 and setting random speeds" <<endl;
					cout << "*** Loading the configuration config.0 ***" << endl << endl;
					M.SetRandomSpeed();
				}
				else if (r1 == 0 && r2 == -1) {
					cerr << "Warning: cannot read the file old.0. Using old.final and setting random speeds instead" <<endl;
					cout << "*** Loading the configuration old.final ***" << endl << endl;
					M.SetRandomSpeed();
				}
				else {
					cout << "*** Loading the configurations old.0 and old.config ***" << endl << endl;
					M.SetSpeedInput();
				}
			}
			else {
				cout << "*** Loading the configuration config.0 ***" << endl << endl;
				r1 = M.SetConfig("config.0");
				if (r1 != 0) {
					cerr << "Error: configuration file config.0 not available. Aborting the simulation" <<endl;
					return -1;
				}
				M.SetRandomSpeed();
			}
			
			// Equilibration cycles
			t1 = M.Equilibrium(NeqCycles, NeqSteps);
			
			// Simulation
			t2 = M.Simulate();
			M.OutputProgMeasure();
			
			cout << "*** Elapsed Times ***" << endl;
			cout << "> Equilibration cycles: " << t1/1000. << " ms" << endl;
			cout << "> Simulation cycles: " << t2/1000. << " ms" << endl<<endl;
			
			// Output of the last configurations
			cout << "> Output (and overwrite, if existing) the final configurations to old.0 and old.final? (y|n) ";
			cin >> option;
			cout << endl << endl;
			if (option == "y" || option == "Y") {
				M.OutputConfig("old.final");
				M.OutputOldConfig("old.0");
			}
		}
		
		cout << "> Start another simulation? (y|n) ";
		cin >> option;
		if (option != "y" && option != "Y") break;
		cout << endl << endl;
	}
	return 0;
}
	
