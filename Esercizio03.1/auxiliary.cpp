#include "auxiliary.h"
#include <cmath>
#include <cstdlib>
#include <chrono>
using namespace std;

double d(double S, double t, double r, double sigma, double K, double T) {
	return (log(S/K)+r+sigma*sigma*(T-t)/2.)/(sigma*sqrt(T-t));
}

double N(double x) {
	return 0.5*(1+erf(x/sqrt(2)));
}

// European Call-option price
double C(double S, double t, double r, double sigma, double K, double T) {
	return S*N(d(S, t, r, sigma, K, T)) - K*exp(-r*(T-t))*N(d(S, t, r, sigma, K, T) -sigma*sqrt(T-t));
}

// European Put-option price
double P(double S, double t, double r, double sigma, double K, double T) {
	return S*(N(d(S, t, r, sigma, K, T))-1) - K*exp(-r*(T-t)) * (N(d(S, t, r, sigma, K, T) - sigma*sqrt(T-t))-1);
}

double OutputData(string filename, vector<double> mean, vector<double> sigma, unsigned int M) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		out << "#Measures " << M << endl;
		out << "#Blocks " << mean.size() << endl;
		for (unsigned int i = 0; i < mean.size(); i++) {
			out << mean[i] << " " << sigma[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
} 

/*
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
*/
