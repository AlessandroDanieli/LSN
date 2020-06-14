#ifndef __Auxiliary__
#define __Auxiliary__
#include <fstream>
#include <string>
#include <vector>
#include "stat.h"
using namespace std;

double OutputProg(string filename, const vector<double>& av_prog, const vector<double>& error_prog, int M);
double OutputData(string filename, Stat s);
void Print(const vector<double> v);

#endif

