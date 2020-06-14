#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include "random.h"
#include "function.h"
#include "stat.h"
#include "auxiliary.h"
using namespace std;

int main(int argc, char** argv) {
	Random rnd;
	rnd.RandGenSetup();
	
	double xmin = 0, xmax = 1;
	unsigned int M = 60000, N = 500;

	Function* f = new Integrand();
	StatIntegral s(&rnd, M, N, f, xmin, xmax); // Data vector, Integral vector
	double t1 = s.FillUniform(); // Uniform sampling
	double t2 = s.BlockIntegrate(); // Calculate each block's integral
	
	// ## Integral and progressive average ## //
	string filename = "prog.out";
	vector<double> average_prog, error_prog;
	Stat prog(s.GetIntegral(), s.GetIntegral().size(), 1); // Data vector, sample_size = size, (one block)
	double t3 = prog.ProgAverage(average_prog, error_prog);
	double ot1 = OutputProg(filename, average_prog, error_prog, M);
	
	// ## Importance sampling - Integral and progressive average ## //
	double m = -1, q = 1; // Linear Importance Sampling parameters
	filename = "prog_IS.out";
	LinearIS* p = new LinearIS(m, q, xmin, xmax); // Linear Importance Sampling p(x) ~ mx+q (* norm)
	double t4 = s.FillLinear(m, q); // Fill data vector with linear distribution p(x) ~ mx+q
	double t5 = s.BlockIntegrate(p);
	
	prog.SetData(s.GetIntegral());
	prog.SetSampleSize(s.GetIntegral().size());
	prog.SetBlocks(1);
	double t6 = prog.ProgAverage(average_prog, error_prog);
	double ot2 = OutputProg(filename, average_prog, error_prog, M);
	
	cout << "Elapsed Times: " << endl;
	cout << "> FillUniform = " << t1 << " \u03BCs" << endl;
	cout << "> BlockIntegrate = " << t2 << " \u03BCs" << endl;
	cout << "> ProgAverage = " << t3 << " \u03BCs" << endl;
	cout << "> FillLinear = " << t4 << " \u03BCs" << endl;
	cout << "> BlockIntegrate(p) = " << t5 << " \u03BCs" << endl;
	cout << "> ProgAverage = " << t6 << " \u03BCs" << endl;
	
	cout << "> OutputProg = " << ot1 << " \u03BCs" << endl;
	cout << "> OutputProg = " << ot2 << " \u03BCs" << endl;
	
	rnd.SaveSeed();
	return 0;
}
