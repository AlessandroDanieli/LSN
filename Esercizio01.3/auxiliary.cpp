#include "auxiliary.h"
using namespace std;

void Print(const vector<double>& v) {
	for (unsigned int i = 0; i < v.size(); i++) {
		cout << v[i] << endl;
	}
	cout << endl;
}

double GenLines(vector<double>& x_lines, double d, double min, double max) {
	auto t1 = chrono::high_resolution_clock::now();
	if (min >= max) {
		cerr << "Error [Auxiliary.cpp :: GenLines]: invalid parameters min, max" << endl;
		exit(-1);
	}
	double n = (max-min)/d;
	if (n - (int)n != 0) {
		cerr << "Error [Auxiliary.cpp :: GenLines]: invalid partition ((max-min)/d is not an integer)" << endl;
		exit(-1);
	}
	x_lines.clear();
	for (unsigned int i = 0; i <= n; i++) {
		x_lines.push_back(min + i*d);
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double GenNeedles(Random* rnd, vector<Needle>& needle_array, unsigned int M, double length, double min, double max) {
	auto t1 = chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < M; i++) {
		needle_array.push_back(Needle(rnd->Uniform(min, max), rnd->Uniform(min, max), length, rnd->Angle()));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double CountNeedles(vector<Needle>& needle_array, vector<double>& x_lines, vector<double>& count, double l, double d) {
	auto t1 = chrono::high_resolution_clock::now();
	count.clear();
	for (unsigned int i = 0; i < needle_array.size(); i++) {
		count.push_back((double)needle_array[i].IsCrossed(x_lines));
	}
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double ConvertToPi(vector<double>& pi_raw, double l, double d) {
	auto t1 = chrono::high_resolution_clock::now();	
	for (unsigned int i = 0; i < pi_raw.size(); i++) {
		if (pi_raw[i] != 0) pi_raw[i] = 2.*l/(d*pi_raw[i]);
		else {
			cerr << "Error [auxiliary.cpp :: ConvertToPi]: measured probability is zero" << endl;
			exit(-1);
		}
	} 
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
}

double OutputData(string filename, vector<double>& prog_average, vector<double>& prog_error, unsigned int M, double l, double d) {
	auto t1 = chrono::high_resolution_clock::now();
	ofstream out(filename);
	if (out.is_open()) {
		out << "#Measures: " << M << endl;
		out << "#Blocks: " << prog_average.size() << endl;
		out << "#NeedleLength: " << l << endl;
		out << "#LineSeparation: " << d << endl;
		for (unsigned int i = 0; i < prog_average.size(); i++) {
			out << i+1 << ' ' << prog_average[i] << ' ' << prog_error[i] << endl;
		}
	}
	out.close();
	auto t2 = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::microseconds>(t2-t1).count();
} 
