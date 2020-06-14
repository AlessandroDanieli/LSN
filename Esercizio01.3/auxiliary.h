#ifndef __Auxiliary__
#define __Auxiliary__
#include "random.h"
#include "needle.h"
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>

void Print(const vector<double>& v);
double GenLines(vector<double>& x_lines, double d, double min, double max);
double GenNeedles(Random* rnd, vector<Needle>& needle_array, unsigned int M, double length, double min, double max);
double CountNeedles(vector<Needle>& needle_array, vector<double>& x_lines, vector<double>& count, double l, double d);
double ConvertToPi(vector<double>& prog_average_raw, double l, double d);
double OutputData(string filename, vector<double>& prog_average, vector<double>& prog_error, unsigned int M, double l, double d);

#endif

