#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include "random.h"
#include "randomwalk.h"
#include "stat.h"
using namespace std;

int main(int argc, char** argv) {
	Random rnd;
	rnd.RandGenSetup();
	
	// RandomWalk setup
	Vector Origin({0, 0, 0});
	Vector x({1, 0, 0}), y({0, 1, 0}), z({0, 0, 1});
	vector<Vector> Generators({x, y, z});
	
	Randomwalk C(&rnd, Origin);
	RandomwalkLattice L(&rnd, Origin, Generators);
	
	// Walk and Output
	unsigned int Nmoves = 100, Nwalks = 1000;
	ofstream Out("walk.out");
	double R_C[Nmoves], R_L[Nmoves], R_C2[Nmoves], R_L2[Nmoves];
	
	if (Out.is_open()) {
		Out << "#Walks " << Nwalks << endl;
		Out << "#Steps " << Nmoves << endl;
		for (unsigned int j = 0; j < Nmoves; j++) {
			R_C[j] = 0;
			R_C2[j] = 0;
			R_L[j] = 0;
			R_L2[j] = 0;
		}
		for (unsigned int i = 0; i < Nwalks; i++) {
			for (unsigned int j = 0; j < Nmoves; j++) {
				C.RandomMove(1);
				R_C[j] += C.Distance(Origin);
				R_C2[j] += C.Distance(Origin)*C.Distance(Origin);
				L.RandomMove();
				R_L[j] += L.Distance(Origin);
				R_L2[j] += L.Distance(Origin)*L.Distance(Origin);
			}
			C.SetPosition(Origin);
			L.SetPosition(Origin);
		}
		for (unsigned int j = 0; j < Nmoves; j++) {
			R_C[j] /= Nwalks;
			R_L[j] /= Nwalks;
			Out << j+1 << " " << R_C[j] << " " << sqrt(R_C2[j]/Nwalks-R_C[j]*R_C[j]) << " " \
				<< R_L[j] << " " << sqrt(R_L2[j]/Nwalks-R_L[j]*R_L[j]) << endl;
		}
	}
	else {
		cerr << "Error: the output file cannot be opened" << endl;
		exit(-1);
	}
	Out.close();
	
	for (int i = 0; i < 100; i++) rnd.Uniform(); // Introduced only in order to produce a better-looking plot
	Nwalks = 4; Nmoves = 100;
	ofstream OutC("continuous.xyz"), OutL("lattice.xyz");
	if (OutC.is_open() && OutL.is_open()) {
		OutC << "#Nwalks " << Nwalks << endl << "Nmoves " << Nmoves << endl;
		OutL << "#Nwalks " << Nwalks << endl << "Nmoves " << Nmoves << endl;
		for (unsigned int n = 0; n < Nwalks; n++) {
			C.SetPosition(Origin); L.SetPosition(Origin);
			OutC << 0 << " " << 0 << " " << 0 << endl;
			OutL<< 0 << " " << 0 << " " << 0 << endl;
			for (unsigned int i = 0; i < Nmoves-1; i++) {
				C.RandomMove(1);
				OutC << C.GetPosition()[0] << " " << C.GetPosition()[1] << " " << C.GetPosition()[2] << endl;
				L.RandomMove();
				OutL<< L.GetPosition()[0] << " " << L.GetPosition()[1] << " " << L.GetPosition()[2] << endl;
			}
			OutC << "#" << endl; 
			OutL << "#" << endl;
		}
		OutC.close(); OutL.close();
	}
	else {
		cerr << "Error: the output file cannot be opened" << endl;
		exit(-1);
	} 
	
	rnd.SaveSeed();
	return 0;
}
