#include "random.h"
#include "simulannealing.h"
#include <iostream>
#include <fstream>
#include <chrono>
//#include "/home/alessandro/Scrivania/LSN/openmpi-4.0.3/ompi/include/mpi.h"
using namespace std;

// Input: parameters of the computation
void Input(const string filename, SimulAnnealing& S, char& c, double& Tmin, double& Tmax, unsigned int& Tcycles) {
	ifstream Input;
	Input.open(filename);
	unsigned int Nx, d;
	double rad;
	Input >> d; 			Input.ignore(120, '\n');	 	
	Input >> Nx; 			Input.ignore(120, '\n');		S.SetParameters(Nx, d);
	Input >> Nx; 			Input.ignore(120, '\n');		S.SetCycles(Nx);
	Input >> d;				Input.ignore(120, '\n'); 		S.SetNorm(d);
	Input >> c; 			Input.ignore(120, '\n'); 		c = tolower(c);
	Input >> rad;			Input.ignore(120, '\n');

	if (c == 'c') {
		S.RandomSphSurface(rad);
	}
	else if (c == 's') {
		S.RandomCube(rad);
	}

	S.GenPath();
	
	string Temp;
	getline(Input, Temp);
	int pos = Temp.find(":");
	int hash_pos = Temp.find("#");
	
	// Case T (one temperature value, one cycle) -> Tcycles set to 0, T = Tmin
	if (pos == -1 || pos > hash_pos) {
		int r = Temp.find(" ");
		while (r != -1 && r < hash_pos) {
			Temp.erase(r, 1);
			r = Temp.find(" ");
		}
		Tmin = stod(Temp.substr(0, hash_pos)); // Single value of T
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
	Input.close();
	c = toupper(c); // Used dynamically in the output filename
	return;
}

// Output: L values of a single temperature cycle
void OutputL(const string filename, vector<double> L, double T) {
	cout << "> Saving L to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#MutationCycles: " << L.size() << endl;
		out << "#Temperature: " << T << endl;
		
		for (unsigned int i = 0; i < L.size(); i++) {
			out << L[i] << endl; 
		}
		out.close();
	}
	else {
		cerr << "Error [OutputL(const string, vector<double>, double)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void OutputParameters(const string filename, SimulAnnealing S, char c, double Tmin, double Tmax, double Tcycles) {
	cout << "> Saving output parameters to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#PointDistribution: ";
		if (c == 'C') out << "Circumference" << endl;
		else if (c == 'S') out << "Square" << endl;
		out << "#Norm: " << S.GetNorm() << endl;
		out << "#TemperatureRange: " << Tmin << " : " << Tmax << endl;
		out << "#TemperatureCycles: " << Tcycles << endl;
		out.close();
	}
	else {
		cerr << "Error [OutputL(const string, vector<double>, double)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}


int main(int argc, char* argv[]) {
	auto t1 = chrono::high_resolution_clock::now();
	
	Random Rnd;
	Rnd.RandGenSetup();	
	
	unsigned int Tcycles;
	double T, Tmin, Tmax;
	char c;
	string option;
	SimulAnnealing S(&Rnd);
	cout << "*** Traveling Salesman Problem ***" << endl;
	cout << "***    Simulated  Annealing    ***" << endl<<endl;
	for (unsigned int iter = 0;  ; iter++) {
		// Read the parameters (point distribution, number of cycles, temperature...)
		Input("input.dat", S, c, Tmin, Tmax, Tcycles); 
		// Output the parameters (referred to a set of future outputs at different T)
		OutputParameters("param"+to_string(S.GetNorm())+string(1, c)+".out", S, c, Tmin, Tmax, Tcycles); 
		// Clean the existing files (referred to the same point distribution and the same norm, which could alter the plots)
		
		vector<double> L(S.GetCycles(), 0);
		
		if (Tcycles != 0) { // Tcycles in [Tmin, Tmax]
			if (Tcycles == 1) { // One cycle at T = mean(Tmin, Tmax) -> Use the following block of code
				Tmin = (Tmin + Tmax)*0.5;
				goto SINGLE_CYCLE; // next block of code: only one cycle at fixed temperature
			}
			else {
				// If Tcycles > 1, the output files are deleted (to avoid undesired "overlaps" of data when plotting)
				int res = system( ("rm -rf PathLengths/L"+to_string(S.GetNorm())+string(1,c)+"*.out").c_str() );
				ofstream out;
				// Decreasing temperature cycles.
				for (unsigned int n = 0; n < Tcycles; n++) {
					T = Tmax + (Tmin-Tmax)*(double)n/(double)(Tcycles-1);
					S.SetT(T);
					L = S.Cycle(0.2); // Returned vector of <Mutation cycles> values of L (one for each generated path)
					OutputL("PathLengths/L"+to_string(S.GetNorm())+string(1, c)+to_string(n+1)+".out", L, T);
				}
			}
		}
		else {
			SINGLE_CYCLE:
			T = Tmin;
			S.SetT(T);
			L = S.Cycle(0.2); // Returned vector of <Mutation cycles> values of L (one for each generated path)
			OutputL("PathLengths/L"+to_string(S.GetNorm())+string(1, c)+to_string(iter+1)+".out", L, T);
		}
		
		S.OutputPathCoord("Path"+to_string(S.GetNorm())+string(1, c)+".out");
	
	
		cout << "> Start another simulation? (y/n) ";
		cin >> option;
		if (option != "y" && option != "Y") break;
		cout << endl << endl;
	}
	Rnd.SaveSeed();
	
	auto t2 = chrono::high_resolution_clock::now();
	cout << "> Elapsed time = " << chrono::duration_cast<chrono::microseconds>(t2-t1).count()/1000000. << " s" << endl;
	return 0;
}
