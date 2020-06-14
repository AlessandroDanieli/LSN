#ifndef __Auxiliary__
#define __Auxiliary__
#include "random.h"
#include "stat.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

double ProgOutput(string filename, vector<double>& prog_average, vector<double>& prog_error, int L);
double ChiOutput(string filename, vector<double>& chi_array, double& chi2, unsigned int M, unsigned int blocks);
double Chi(const Stat& s, vector<double>& chi_array, double& chi2, unsigned int M, unsigned int blocks);
#endif

