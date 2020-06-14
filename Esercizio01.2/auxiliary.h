#ifndef __Auxiliary__
#define __Auxiliary__
#include "random.h"
#include "stat.h"
#include <iostream>
#include <fstream>
#include <string>

double RollDice(Random rnd, Stat& data_unif, Stat& data_exp, Stat& data_lorentz, int M, \
	double lambda, double mean_lorentz, double gamma_lorentz);
double RollDiceContinuous(Random* , Stat& , Stat& , Stat& , unsigned int , double , double , double );

double Analyse(Stat& data_unif, Stat& data_exp, Stat& data_lorentz, \
	vector<double>* dice_std, vector<double>* dice_exp, vector<double>* dice_lorentz, unsigned int n_set, unsigned int* L, unsigned int N);

double OutputData(string filename, Stat& data_unif, Stat& data_exp, Stat& data_lorentz, unsigned int M);
double OutputAnalysis(string filename, vector<double>* dice_std, vector<double>* dice_exp, vector<double>* dice_lorentz, unsigned int M, unsigned int N, unsigned int n_set, unsigned int* L);

#endif

