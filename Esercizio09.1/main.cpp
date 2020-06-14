#include "random.h"
#include "genetic.h"
#include <iostream>
#include <chrono>
#include <cmath>
using namespace std;

void Input(const string filename, Genetic& G, char& c, unsigned int& Ngen) {
	ifstream Input;
	Input.open(filename);
	unsigned int Npop, Nx, d;
	double r;
	Input >> d; 			Input.ignore(120, '\n');	 	
	Input >> Nx; 			Input.ignore(120, '\n');
	Input >> Npop; 			Input.ignore(120, '\n');
	G.SetParameters(Npop, Nx, d);
	Input >> Ngen;			Input.ignore(120, '\n');	
	Input >> c; 			Input.ignore(120, '\n'); 		c = tolower(c);
	Input >> r;				Input.ignore(120, '\n');
	
	if (c == 'c') {
		G.RandomSphSurface(r);
	}
	else if (c == 's') {
		G.RandomCube(r);
	}
	G.GenPopulation();
	Input >> d;				Input.ignore(120, '\n'); 		G.SetNorm(d);
	c = toupper(c);
	Input.close();
	return;
}

void OutputL(const string filename, vector<double> L, vector<double> MeanL) {
	cout << "> Saving <L> to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#Generations: " << MeanL.size()<< endl;
		out << "#Individuals: " << L.size()<< endl;
		for (unsigned int i = 0; i < MeanL.size(); i++) {
			out << MeanL[i] << endl; 
		}
		out << "#" << endl;
		for (unsigned int i = 0; i < L.size()-1; i++) {
			out << L[i] << endl; 
		}
		out << L[L.size()-1];
		out.close();
	}
	else {
		cerr << "Error [OutputMeanL(const string, vector<double>, vector<double>)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

int main() {
	auto t1 = chrono::high_resolution_clock::now();
	
	Random Rnd;
	Rnd.RandGenSetup();	
	Genetic G;
	G.SetRandGen(&Rnd);
	unsigned int Ngen;
	char c;
	Input("input.dat", G, c, Ngen);
	vector<double> L, MeanL;
	
	for (int cycle = 0; cycle < Ngen; cycle++) {
		cout << "*** Generation " << cycle+1 << " ***" << endl;
		G.Measure();
		G.Sort();
		G.Select();
		G.PrintBestL();
		MeanL.push_back(G.MeanL());
	}
	L = G.GetL(); // Best L as a function of the population
	cout << "> Best path length = " << G.GetBestL() << endl;
	G.SetBestPath();
	G.OutputBestPath("path"+string(1, c)+to_string(G.GetNorm())+".out");
	OutputL("L"+string(1, c)+to_string(G.GetNorm())+".out", L, MeanL);
	Rnd.SaveSeed();
	
	auto t2 = chrono::high_resolution_clock::now();
	cout << "> Elapsed time = " << chrono::duration_cast<chrono::microseconds>(t2-t1).count()/1000000. << " s" << endl;
	return 0;
}
