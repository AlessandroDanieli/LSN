#ifndef __Auxiliary__
#define __Auxiliary__
#include "random.h"
#include "stat.h"
#include <iostream>
#include <fstream>
#include <string>

double d(double S, double t, double r, double sigma, double K, double T);
double N(double x);

// European Call-option price
double C(double S, double t, double r, double sigma, double K, double T);
// European Put-option price
double P(double S, double t, double r, double sigma, double K, double T);

double OutputData(string filename, vector<double> mean, vector<double> sigma, unsigned int M);


#endif

