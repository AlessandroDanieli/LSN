#include "auxiliary.h"
#include <cmath>
#include <cstdlib>
#include <chrono>
using namespace std;

// Continuous, Uniform(0, 6), Exp and Lorentz unbound
double RollDiceContinuous(Random* rnd, Stat& data_unif, Stat& data_exp, Stat& data_lorentz, unsigned int M, \
	double lambda, double mean_lorentz, double gamma_lorentz) {
	auto t1 = chrono::high_resolution_clock::now();
	
	for (unsigned int i = 0; i < M; i++) {
		data_unif.Append(rnd->Uniform(0, 6));
		data_exp.Append(rnd->Exp(lambda));
		data_lorentz.Append(rnd->Lorentz(mean_lorentz, gamma_lorentz));
	}
	
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

// Selects only the integers between 1 and 6
double RollDice(Random rnd, Stat& data_unif, Stat& data_exp, Stat& data_lorentz, unsigned int M, \
	double lambda, double mean_lorentz, double gamma_lorentz) {
	auto t1 = chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < M; i++) {
		data_unif.Append(ceil(rnd.Uniform()*6.));
	}
	
	double p;
	unsigned int k = 0;
	while (k < M) {
		p = rnd.Exp(lambda);
		if (p >= 0 && p < 6) {
			data_exp.Append(ceil(p));
			k++;
		}	
	}
	k = 0;
	while (k < M) {
		p = rnd.Lorentz(mean_lorentz, gamma_lorentz);
		if (p >= 0 && p < 6) {
			data_lorentz.Append(ceil(p));
			k++;
		}	
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Analyse(Stat& data_unif, Stat& data_exp, Stat& data_lorentz, \
vector<double>* dice_std, vector<double>* dice_exp, vector<double>* dice_lorentz, unsigned int n_set, unsigned int* L, unsigned int N) {
	auto t1 = chrono::high_resolution_clock::now();
	data_unif.SetBlocks(N);
	data_exp.SetBlocks(N);
	data_lorentz.SetBlocks(N);
	for (unsigned int i = 0; i < n_set; i++) {
		if (L[i]*N > data_unif.GetSampleSize()) {
			cerr << "Error [Analyse :: auxiliary.cpp]: <data> index out of range" << endl;
			exit(-1);
		}
	}
	
	for (unsigned int i = 0; i < n_set; i++) {
		data_unif.SetSampleSize(N*L[i]); // Sets the size of the data set for a given block number and length.
		data_exp.SetSampleSize(N*L[i]);
		data_lorentz.SetSampleSize(N*L[i]);
		data_unif.BlockAverage(dice_std[i]);
		data_exp.BlockAverage(dice_exp[i]);
		data_lorentz.BlockAverage(dice_lorentz[i]);
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double OutputData(string filename, Stat& data_unif, Stat& data_exp, Stat& data_lorentz, unsigned int M) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		for (unsigned int i = 0; i < M; i++) {
			out << data_unif[i] << ' ' << data_exp[i] << ' ' << data_lorentz[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
} 

double OutputAnalysis(string filename, vector<double>* dice_std, vector<double>* dice_exp, vector<double>* dice_lorentz, \
	unsigned int M, unsigned int N, unsigned int n_set, unsigned int* L) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		//Output Parameters
		out << "#Measures: " << M << endl;
		out << "#Blocks: " << N << endl;
		out << "#Sets: " << n_set << endl;
		out << "#Block_lengths: ";
		for (unsigned int i = 0; i < n_set-1; i++) {
			out << L[i] << ' ';
		}
		out << L[n_set-1] << endl;
	
		// Output dice measures
		for (unsigned int i = 0; i < n_set; i++) {
			for (unsigned int j = 0; j < N; j++) {
				out << dice_std[i][j] << ' ' << dice_exp[i][j] << ' ' << dice_lorentz[i][j] << endl;
			}
			out << endl; // Separation between sets with "\n"
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
