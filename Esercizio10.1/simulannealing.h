#ifndef __TSP__
#define __TSP__
#include "random.h"
#include "vector.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

class SimulAnnealing {
private:
	Random *Rnd;
	vector<Vector> X; // vector of "Vector" objects (d-dimensional points)
	unsigned int Nx, d, Norm, Cycles; // Number of points, dimension of the space
	vector<unsigned int> Path;
	double L, T;
	vector<unsigned int> FinalPath; // Best path to be printed on file
	
public:
	// Constructors/Destructor
	SimulAnnealing();
	SimulAnnealing(Random* r);
	SimulAnnealing(int _N, int _d);
	SimulAnnealing(Random* r, int _N, int _d);
	~SimulAnnealing(); 
	
	// Get Methods
	unsigned int GetNx() const { return Nx; }
	unsigned int GetD() const { return d; }
	unsigned int GetNorm() const { return Norm; }
	unsigned int GetCycles() const { return Cycles; }
	double GetL() const { return L; }
	Random* GetRandGen() { return Rnd; }
	
	// Set Methods
	void SetParameters(int _Nx, int _d);
	void SetNorm(unsigned int p);
	void SetT(double _T);
	void SetRandGen(Random* r) { Rnd = r; }
	void SetCycles(int N);
	
	// Measure Algorithms - Path
	void Measure();
	void L1();
	void L2();
	
	// Measure Algorithms - NewPath
	double Measure(vector<unsigned int>& NewPath);
	double L1(vector<unsigned int>& NewPath);
	double L2(vector<unsigned int>& NewPath);
		
	// Individuals' Generation || Points' Pattern Generation
	void GenPath();
	void RandomSphSurface(double R);
	void RandomCube(double l);
	
	// Mutation Algorithms
	vector<double> Cycle(double threshold);
	void Mutate(double threshold);
	int Swap(vector<unsigned int>& NewPath, unsigned int i, unsigned int j); 
	int Permute(vector<unsigned int>& NewPath, int n_block, int i, int j=-1); // Block permutation (generic)
	int Shift(vector<unsigned int>& NewPath, int i, int n_shift, int n_places);
	int Invert(vector<unsigned int>& NewPath, unsigned int i, unsigned int n_inv);
	//void CrossingOver(unsigned int p, unsigned int q, unsigned int cut);
	
	// Screen Output Methods
	void PrintPath() const;
	void PrintPath(vector<unsigned int> NewPath) const; 
	void PrintX() const;
	/*
	void PrintIndividual(unsigned int i) const;
	void PrintL() const;
	void PrintOrder() const;
	void PrintBest() const;
	void PrintBestL() const;
	void PrintV(vector<unsigned int> V) const;
	*/
	
	// File Input/Output
	void LoadPoints(const string filename);
	void LoadPathIndices(const string filename);
	void OutputPoints(const string filename) const;
	void OutputPathCoord(const string filename) const;
	void OutputPathIndices(const string filename) const;
	
};

#endif
