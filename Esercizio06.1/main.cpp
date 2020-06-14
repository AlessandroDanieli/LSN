#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include "random.h"
#include "ising.h"
#include "stat.h"
using namespace std;

bool FileExists(const string& filename) {
	struct stat buffer;   
	return (stat(filename.c_str(), &buffer) == 0); 
}

// Start simulation cycles from Tmax to Tmin, using as starting configuration the previous (less equilibration cycles)
void Input(const string filename, Ising& I, bool& ReadConfig, bool* OutputOption, unsigned int Nobs, double& Tmin, double& Tmax, int& Tcycles) {
	ifstream Input;
	Input.open(filename);
	unsigned int t, t2;
	double d;
	char c;
	if (Input.is_open()) {
		// PHYSICAL PARAMETERS
		Input.ignore(250, '\n');	Input.ignore(250, '\n');		Input.ignore(250, '\n');
		
		// Number of spins
		Input >> t; 				Input.ignore(250, '\n');	 	I.SetNspin(t);
		
		// TEMPERATURE
		// In order to read the input line "Tmin:Tmax:Tcycles #[Text...]", the characters ":" and "#" must be located
		string Temp;
		getline(Input, Temp);
		int pos = Temp.find(":");
		int hash_pos = Temp.find("#");
		
		// Case T: one temperature value, one cycle
		if (pos == -1 || pos > hash_pos) {
			int r = Temp.find(" ");
			while (r != -1 && r < hash_pos) {
				Temp.erase(r, 1);
				r = Temp.find(" ");
			}
			Tmin = stod(Temp.substr(0, hash_pos));
			Tcycles = -1; // Used to denote "one single simulation"
		}
		
		// Case Tmin:Tmax:Cycles -> repeating the simulation <Cycles> times between Tmin and Tmax
		else {
			int r = Temp.find(" ");
			while (r != -1 && r < hash_pos) {
				Temp.erase(r, 1);
				r = Temp.find(" ");
			}
			Tmin = stod(Temp.substr(0, pos));
			Temp.erase(0, pos+1);
			pos = Temp.find(":");
			Tmax = stod(Temp.substr(0, pos));
			Temp.erase(0, pos+1);
			pos = hash_pos;
			Tcycles = stoi(Temp.substr(0, pos));
			
			if (Tmin >= Tmax) {
				cerr << "Error [Main :: Input]: Invalid temperature range" << endl;
				exit(-1);
			}
			if (Tcycles <= 0) {
				cerr << "Error [Main :: Input]: Invalid number of cycles" << endl;
				exit(-1);
			}
		}
		// J, external H
		Input >> d; 					Input.ignore(250, '\n');		I.SetJ(d);
		Input >> d; 					Input.ignore(250, '\n');		I.SetHext(d);
		
		// SIMULATION PARAMETERS (Number of steps and blocks, Algorithm option)
		Input.ignore(250, '\n');		Input.ignore(250, '\n');		Input.ignore(250, '\n');
		// Set Nblocks, Nsteps (Stat object)
		Input >> t; 					Input.ignore(250, '\n'); 		
		Input >> t2; 					Input.ignore(250, '\n'); 		I.SetStat(t, t2);
		// Set algorithm option 
		Input >> c; 					Input.ignore(250, '\n'); 		c = tolower(c); 		I.SetAlgOption(c);
		
		// INPUT/OUTPUT
		Input.ignore(250, '\n');		Input.ignore(250, '\n');		Input.ignore(250, '\n');
		Input >> ReadConfig;			Input.ignore(250, '\n');
		for (unsigned int i = 0; i < Nobs; i++) { // Four observables
			Input >> OutputOption[i];		Input.ignore(250, '\n'); 
		}
		
		Input.close();
 	}
 	else {
		cerr << "Error [Main :: Input]: cannot open the file <" << filename << ">" << endl;
		exit(-1);
	}
	return;
}

int main() {
	unsigned int Nobs = 4; // Number of observables: U, C, M, X
	bool OutputOption[Nobs], ReadConfig; // Output observables, read starting configuration
	string option;
	cout << "*** Classic 1D Ising model - Nearest neighbour interaction ***" << endl;
	cout << "***                 Using k_B = 1, \u03BC_B = 1                 ***" << endl << endl;

	Random Rnd;
	Rnd.RandGenSetup();
	Ising I(&Rnd); // Ising object. 4 observables U, C, M, X (per spin) by default

	double Tmin, Tmax;
	int Tcycles;
	
	// Simulation cycle
	for (unsigned int SimCycle = 0;  ; SimCycle++) {
		// Read system parameters and output options
		Input("input.dat", I, ReadConfig, OutputOption, Nobs, Tmin, Tmax, Tcycles);
		
		// Output option (overwrite / append / clean) for each observable's output file
		cout << "> Overwrite the existing output files (otherwise append)? (y/n) ";
		cin >> option;
		cout << endl;
		if (option == "y" || option == "Y") {
			for (unsigned int obs = 0; obs < Nobs; obs++) if (OutputOption[obs]) I.CleanOutput(obs);
		}
		
		// Initial configuration
		if (ReadConfig) I.LoadConfig("config.final"); 	// Read an old configuration
		else I.SetRandomSpin(); 						// Random initial spin configuration
		// Simulation
		cout <<"[input.dat] "<< Tcycles <<" cycles of simulation in the temperature range ["<<Tmin<< ", "<<Tmax<<"]"<<endl<<endl;
		//I.PrintParameters();
		
		// Temperature cycles
		if (Tcycles != -1) { 					// Tcycles in [Tmin, Tmax]
			if (Tcycles == 1) { 				// One cycle at T = mean(Tmin, Tmax) -> Use the next block of code
				Tmin = (Tmin + Tmax)*0.5;
				goto DEFAULT_CYCLE; 			// next block of code [label]
			}
			for (unsigned int cycle = 0; cycle < Tcycles; cycle++) { // Starting from Tmax -> Tmin in <Tcycle> cycles
				I.SetT(Tmax + (Tmin-Tmax)*(double)cycle/(double)(Tcycles-1));
				cout << "> Cycle " << cycle+1 << " [T = " << I.GetT() << "]" << endl;
				
				I.AutoEquilibration(50, 0.01); // Automatic equilibration using the variance of the mean energy U and a threshold
				I.Simulate();
				for (unsigned int obs = 0; obs < Nobs; obs++) {
					if (OutputOption[obs]) {
						I.ComputeObs(obs);
						I.OutputMeasure(obs);
					}
				}
				
				// Reset the internal data storage (each step's U and total Spin)
				I.ClearMeasures(); 
			}
		}
		
		// Single value of T - Single simulation
		else { 					// One temperature value T, one simulation
			DEFAULT_CYCLE: 		// GOTO label, used if Tcycles = 1 (special case) with temperature T = (Tmin+Tmax)/2
				I.SetT(Tmin);
				
				I.AutoEquilibration(50, 0.01);
				I.Simulate();
				for (unsigned int obs = 0; obs < Nobs; obs++) {
					if (OutputOption[obs]) {
						I.ComputeObs(obs);
						I.OutputMeasure(obs);
					}
				}
				I.ClearMeasures(); // Reset the internal data storage (each step's U and total Spin)
		}
		
		// Final configuration - Output
		cout << "> Save the actual spin configuration on <config.final>? (y|n) ";
		cin >> option;
		cout << endl;
		if (option == "y" || option == "Y") I.OutputConfig("config.final");
		
		cout << "> Start another simulation? (y/n) ";
		cin >> option;
		if (option != "y" && option != "Y") break;
		cout << endl << endl;
	}
	
	Rnd.SaveSeed();
	return 0;
}


/*

// Start simulation cycles from Tmax to Tmin, using as starting configuration the previous (less equilibration cycles)
void Input(const string filename, Ising& I, bool* OutputOption, unsigned int Nobs, double& Tmin, double& Tmax, int& Tcycles) {
	ifstream Input;
	Input.open(filename);
	unsigned int t;
	double d;
	char c;
	if (Input.is_open()) {
		Input >> t; 			Input.ignore(250, '\n');	 	I.SetNspin(t);
		Input >> t; 			Input.ignore(250, '\n'); 		I.SetNblocks(t);
		Input >> t; 			Input.ignore(250, '\n'); 		I.SetNsteps(t);
		
		// In order to read the input line "Tmin:Tmax:Tcycles #[Text...]", the characters ":" and "#" must be located
		string Temp;
		getline(Input, Temp);
		int pos = Temp.find(":");
		int hash_pos = Temp.find("#");
		
		// Case T (one temperature value, one cycle)
		if (pos == -1 || pos > hash_pos) {
			int r = Temp.find(" ");
			while (r != -1 && r < hash_pos) {
				Temp.erase(r, 1);
				r = Temp.find(" ");
			}
			Tmin = stod(Temp.substr(0, hash_pos));
			Tcycles = 0;
		}
		
		// Case Tmin:Tmax:Cycles -> repeating the simulation <Cycles> times between Tmin and Tmax
		else {
			int r = Temp.find(" ");
			while (r != -1 && r < hash_pos) {
				Temp.erase(r, 1);
				r = Temp.find(" ");
			}
			Tmin = stod(Temp.substr(0, pos));
			Temp.erase(0, pos+1);
			pos = Temp.find(":");
			Tmax = stod(Temp.substr(0, pos));
			Temp.erase(0, pos+1);
			pos = hash_pos;
			Tcycles = stoi(Temp.substr(0, pos));
			
			if (Tmin >= Tmax) {
				cerr << "Error [Main :: Input]: Invalid temperature range" << endl;
				exit(-1);
			}
			if (Tcycles <= 0) {
				cerr << "Error [Main :: Input]: Invalid number of cycles" << endl;
				exit(-1);
			}
		}
	
		Input >> d; 					Input.ignore(250, '\n');	I.SetJ(d);
		Input >> d; 					Input.ignore(250, '\n');	I.SetHext(d);
		Input >> c; 					Input.ignore(250, '\n'); 	c = tolower(c); 	I.SetAlgOption(c);
		
		// Output Option for each observable
		for (unsigned int i = 0; i < Nobs; i++) { 
			Input >> OutputOption[i];		Input.ignore(250, '\n'); 
		}
		Input.close();
 	}
 	else {
		cerr << "Error [Main :: Input]: cannot open the file <" << filename << ">" << endl;
		exit(-1);
	}
	return;
}

int main() {
	unsigned int Nobs = 4; // Number of observables: U, C, M, X
	bool OutputOption[Nobs];
	string option;
	cout << "*** Classic 1D Ising model - Nearest neighbour interaction ***" << endl;
	cout << "***                 Using k_B = 1, \u03BC_B = 1                 ***" << endl << endl;

	Random Rnd;
	Rnd.RandGenSetup();
	Ising I(&Rnd); // Ising object. 4 observables U, C, M, X (per spin) by default

	double Tmin, Tmax;
	int Tcycles;
	
	// Simulation cycle
	for (unsigned int cycle = 0;  ; cycle++) {
		// Read system parameters and output options
		Input("input.dat", I, OutputOption, Nobs, Tmin, Tmax, Tcycles);
		
		// Output option (overwrite / append / clean) for each observable's output file
		cout << "> Overwrite the existing output files? (y/n) ";
		cin >> option;
		if (option == "y" || option == "Y") 
			for (unsigned int obs = 0; obs < Nobs; obs++) if (OutputOption[obs]) I.CleanOutput(obs);
		cout << endl;
		
		cout <<"[input.dat] "<<Tcycles<< " cycles of simulation in the temperature range ["<<Tmin<< ", "<<Tmax<<"]"<<endl<<endl;
		I.PrintParameters();
		// Simulation
		if (Tcycles != 0) { // Tcycles in [Tmin, Tmax]
			if (Tcycles == 1) { // One cycle at T = mean(Tmin, Tmax) -> Use the following block of code
				Tmin = (Tmin + Tmax)*0.5;
				goto DefaultCycle; // next block of code
			}
			for (unsigned int k = 0; k < Tcycles; k++) { // Starting from Tmax -> Tmin in <Tcycle> cycles
				I.SetT(Tmax + (Tmin-Tmax)*(double)k/(double)(Tcycles-1));
				cout << "> Cycle " << k+1 << " [Temperature = " << I.GetT() << "]" << endl;
				if (k == 0) I.Initialize();
				else I.LoadConfig("config.final"); // Loading the last configuration (sparing some Equilibration cycles)
				I.AutoEquilibration(50, 0.01); // Automatic equilibration using the variance of the mean energy U and a threshold
				
				I.Simulate();
				for (unsigned int obs = 0; obs < Nobs; obs++) {
					if (OutputOption[obs]) {
						I.ComputeObs(obs);
						I.OutputMeasure(obs);
					}
				}
				
				I.OutputConfig("config.final");
				I.ClearMeasures(); // Reset the internal data storage (each step's U and total Spin)
				cout << endl;
			}
		}
		else { // One value of temperature T, one cycle of simulation
			DefaultCycle: // GOTO label, used if Tcycles = 1 (special case) with temperature T = (Tmin+Tmax)/2
				I.SetT(Tmin);
				I.Initialize();
				I.AutoEquilibration(50, 0.01);
				I.Simulate();
				for (unsigned int obs = 0; obs < Nobs; obs++) {
					if (OutputOption[obs]) {
						I.ComputeObs(obs);
						I.OutputMeasure(obs);
					}
				}
				I.ClearMeasures(); // Reset the internal data storage (each step's U and total Spin)
		}

		cout << "> Start another simulation? (y/n) ";
		cin >> option;
		if (option != "y" && option != "Y") break;
		cout << endl << endl;
	}
	
	Rnd.SaveSeed();
	return 0;
}
*/
