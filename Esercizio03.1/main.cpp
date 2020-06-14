#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "brownian.h"
#include "stat.h"
#include "auxiliary.h" // -> Output functions, auxiliary functions
using namespace std;

void print(vector<double> x) {
	for (int i = 0; i < x.size(); i++) cout << x[i] << endl;
}

int main(int argc, char** argv) {
	auto t1 = chrono::high_resolution_clock::now();	
	Random rnd;
	rnd.RandGenSetup();
	
	//unsigned int M = 1000, Nblocks = 200, Nsteps = 100;
	unsigned int M = 100000, Nblocks = 200, Nsteps = 100;
	double r = 0.1, sigma = 0.25;
	double S0 = 100, T = 1, K = 100;
	
	GBM S_direct(&rnd, r, sigma, S0); // Direct sampling of S(T)
	GBM S_discrete(&rnd, r, sigma, S0); // Discrete sampling of S(T), using a specific recursive formula
	rnd.SaveSeed();
	
	Stat C_dis(M, Nblocks), C_dir(M, Nblocks), P_dis(M, Nblocks), P_dir(M, Nblocks);
	
	for (unsigned int m = 0; m < M; m++) {
		S_direct.Set(0, S0);
		S_direct.Evolve(T);
		S_discrete.Set(0, S0);
		S_discrete.Evolve((double)T/Nsteps, Nsteps);
		
		C_dir.Append(C(S_direct.GetS(), T, r, sigma, K, T));
		P_dir.Append(P(S_direct.GetS(), T, r, sigma, K, T));
		C_dis.Append(C(S_discrete.GetS(), T, r, sigma, K, T));
		P_dis.Append(P(S_discrete.GetS(), T, r, sigma, K, T));
	}
	
	vector <double> C_dis_mean, C_dis_sigma, C_dir_mean, C_dir_sigma, P_dis_mean, P_dis_sigma, P_dir_mean, P_dir_sigma;
	C_dir.BlockProgAverageSigma(C_dir_mean, C_dir_sigma);
	C_dis.BlockProgAverageSigma(C_dis_mean, C_dis_sigma);
	P_dir.BlockProgAverageSigma(P_dir_mean, P_dir_sigma);
	P_dis.BlockProgAverageSigma(P_dis_mean, P_dis_sigma);
	
	OutputData("call_direct.out", C_dir_mean, C_dir_sigma, M);
	OutputData("call_discrete.out", C_dis_mean, C_dis_sigma, M);
	OutputData("put_direct.out", P_dir_mean, P_dir_sigma, M);
	OutputData("put_discrete.out", P_dis_mean, P_dis_sigma, M);
	

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Elapsed time = " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;
	return 0;
}
