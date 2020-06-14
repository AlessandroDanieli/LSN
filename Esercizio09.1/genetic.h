#ifndef __Genetic__
#define __Genetic__
#include "random.h"
#include "vector.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

class Genetic {
private:
	Random *Rnd;
	vector<Vector> X; // vector of "Vector" objects (d-dimensional points)
	vector<vector<unsigned int>> Population; // Array of Npop Individuals [= vector of Nx indices]
	vector<double> L; // Weight/cost of each individual
	vector<unsigned int> Order; // Stores the ordered sequence of individuals' indices (avoiding expensive exchanges)
	vector<unsigned int> FinalPath; // Best path to be printed on file
	unsigned int Npop, Nx, d, Norm; // Number of points, dimension of the space
	
public:
	// Constructors/Destructor
	Genetic();
	Genetic(int _Npop, int _N, int _d);
	Genetic(Random* r, int _Npop, int _N, int _d);
	~Genetic(); 
	
	// Get Methods
	unsigned int GetNpop() const { return Npop; }
	unsigned int GetNx() const { return Nx; }
	unsigned int GetD() const { return d; }
	unsigned int GetNorm() const { return Norm; }
	double GetBestL() { return L[0]; }
	vector<double> GetL() { return L; }
	Random* GetRandGen() { return Rnd; }
	
	// Set Methods
	void SetParameters(int _Npop, int _Nx, int _d);
	void SetNorm(unsigned int p);
	void SetRandGen(Random* r) { Rnd = r; }
	void SetBestPath(); // Stores the best path in vector<unsigned int> FinalPath
	//void SetNpop(int _Npop);
	//void SetNx(int _Nx);
	//void SetD(int _d);
	
	// Measure Algorithms
	void Measure();
	void L1();
	void L2();
	double MeanL() const;
	
	// Sort the vectors L, Order, keeping FIXED the position of the individuals in Population[][]
	void Sort(int sx = -1, int dx = -1);
	
	// Individuals' Generation || Points' Pattern Generation
	void GenPopulation();
	void RandomSphSurface(double R);
	void RandomCube(double l);
	
	// Mutation Algorithms
	void Mutate(int first = -1, int last = -1);
	int Swap(unsigned int p, unsigned int i, unsigned int j); 
	int Permute(unsigned int p, int n_block, int i, int j=-1); // Block permutation (generic)
	int Shift(unsigned int p, int i, int n_shift, int n_places);
	int Invert(unsigned int p, unsigned int i, unsigned int n_inv);
	void CrossingOver(unsigned int p, unsigned int q, unsigned int cut);
	
	// Auxiliary algorithms ( Used in CrossingOver() )
	int BinarySearch(vector<unsigned int> V, unsigned int first, unsigned int last, unsigned int find) const;
	void Shift(vector<unsigned int>& V, int i, int i_new) const;
	
	// Selection Algorithm
	void Select();
	
	// Screen Output Methods
	void PrintX() const;
	void PrintIndividual(unsigned int i) const;
	void PrintPopulation() const;
	void PrintL() const;
	void PrintOrder() const;
	void PrintBest() const;
	void PrintBestL() const;
	void PrintV(vector<unsigned int> V) const;
	
	// File Input/Output
	void LoadPoints(const string filename);
	void OutputPoints(const string filename) const;
	void OutputBestPath(const string filename) const;
};

#endif
