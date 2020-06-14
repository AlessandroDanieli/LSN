#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "random.h"
#include "stat.h"
#include "metropolisorbitals.h"
using namespace std;

bool FileExists(const string& filename);
int Input(string filename, int**& QN, unsigned int& Norbitals, unsigned int& Nsample, unsigned int& Nblocks, Vector*& Start, double*& RangeOfMotion, bool*& CalibOption, unsigned int& CalibSteps, unsigned int Nsampling, unsigned int& Sample_first, unsigned int& Sample_last);
//int InputQN(const string filename, int**& QN, unsigned int& Norbitals);

void OutputFinalSetup(string file, int** QN, unsigned int Norb, unsigned int Ns, unsigned int Nbl, Vector* St, double* R, bool* CalibOption, unsigned int CalibSteps, unsigned int Sample_first, unsigned int Sample_last);
void OutputProgRadius(vector<double> prog_average, vector<double> prog_sigma, string filename, unsigned int Nsample);

void PrintQN(int**& QN, unsigned int& Norbitals);
void PrintRange(double* RangeOfMotion, unsigned int Norbitals);
