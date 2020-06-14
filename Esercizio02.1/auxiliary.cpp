#include "auxiliary.h"
#include <iostream>
using namespace std;

double OutputData(string filename, Stat s) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		for (unsigned int i = 0; i < s.GetDataSize(); i++) {
			out << s[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
} 

double OutputProg(string filename, const vector<double>& av_prog, const vector<double>& error_prog, int M) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		out << "#Measures: " << M << endl;
		out << "#Blocks: " << av_prog.size() << endl;
		for (unsigned int i = 0; i < av_prog.size(); i++) {
			out << i+1 << ' ' << av_prog[i] << ' ' << error_prog[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

void Print(const vector<double> v) {
	for (unsigned int i = 0; i < v.size(); i++) {
		cout << v[i] << endl;
	}
	cout << endl;
}
