#include "auxiliary.h"
#include <cmath>
#include <chrono>
using namespace std;

double ProgOutput(string filename, vector<double>& prog_average, vector<double>& prog_error, int L) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		for(unsigned int i = 0; i < prog_average.size(); i++) {
		  out << (i+1)*L << " " << prog_average[i] << " " << prog_error[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double Chi(const Stat& s, vector<double>& chi_array, double& chi2, unsigned int M, unsigned int L) {
	auto t1 = chrono::high_resolution_clock::now();
	if (M*L > s.GetDataSize()) {
		cerr << "Error [auxiliary.cpp :: Chi(const Stat&, vector<double>&, unsigned int, unsigned int]: chi sample size out of range" << endl;
		exit(-1);
	}
	//unsigned int L = (unsigned int) M/blocks;
	//double exp_value = (double)L;
	vector<double> count (M, 0);
	
	chi_array.clear();
	for (unsigned int i = 0; i < M*L; i++) {
		count[(int) ( s[i] * M ) ]++; // Valid since the values belong to (0, 1)
	}
	chi2 = 0;
	double exp = (double)L;
	for (unsigned int k = 0; k < M; k++) {
		chi_array.push_back((count[k] - exp)*(count[k] - exp)/exp);
		chi2 += chi_array[k];
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double ChiOutput(string filename, vector<double>& chi_array, double& chi2, unsigned int M, unsigned int L) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		out << M << " " << L << " " << chi2 << endl;
		for(unsigned int i = 0; i < M; i++){
		  out << chi_array[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

/*
double Chi(const Stat& s, vector<double>& chi_array, unsigned int M, unsigned int blocks) {
	auto t1 = chrono::high_resolution_clock::now();
	if (M > s.GetDataSize()) {
		cerr << "Error [auxiliary.cpp :: Chi(const Stat&, vector<double>&, unsigned int, unsigned int]: chi sample size out of range" << endl;
		exit(-1);
	}
	if (M % blocks != 0) {
		cerr << "Error [auxiliary.cpp :: Chi(const Stat&, vector<double>&, unsigned int, unsigned int]: M % blocks != 0" << endl;
		exit(-1);
	}
	unsigned int L = (unsigned int) M/blocks;
	double exp_value = (double)L/blocks;
	vector<double> count (blocks, 0);
	
	chi_array.clear();
	for (unsigned int i = 0; i < blocks; i++) {
		chi_array.push_back(0);
		for (unsigned int j = 0; j < blocks; j++) {
			count[j] = 0;
		}
		// Sum over each data block. count is referred to a single data block (i index)
		for (unsigned int j = 0; j < L; j++) {
			count[(int) ( s[i*L+j] * blocks ) ]++; // Valid since the values belong to (0, 1)
		}
		// New chi_array value for the i-th block
		for (unsigned int k = 0; k < blocks; k++) {
			chi_array[i] += (count[k] - exp_value)*(count[k] - exp_value)/exp_value;
		}
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
double ChiOutput(string filename, vector<double>& chi_array, unsigned int M, unsigned int blocks) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		out << M << ' ' << blocks << endl;
		for(unsigned int i = 0; i < blocks; i++){
		  out << chi_array[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}
*/

