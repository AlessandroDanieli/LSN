#include "random.h"
#include "stat.h"
#include "auxiliary.h"
#include "position.h"
#include "needle.h"
#include <chrono>

using namespace std;

int main (int argc, char *argv[]) {
	Random rnd;
	rnd.RandGenSetup();

	unsigned int M = 1000000, N = 500; // Number of needles, number of blocks (pi estimations)
	double min = -100, max = 100, l = 12, d = 25;
	vector<double> x_lines;
	GenLines(x_lines, d, min, max);
	
	vector<Needle> needle_array;
	double t1 = GenNeedles(&rnd, needle_array, M, l, min, max);
	
	vector<double> count, prog_average, prog_error;;
	double t2 = CountNeedles(needle_array, x_lines, count, l, d);
	Stat pi(count, M, N);
	double t3 = pi.BlockAverage(prog_average); // Average probability for each block
	double t4 = ConvertToPi(prog_average, l, d); // pi = 2*L/(P*d)
	double t5 = pi.CopyData(prog_average); // Reusing Stat pi
	double t6 = pi.BlockProgAverage(prog_average, prog_error); // Progressive average, error
	
	string filename = "buffon.out";
	double ot1 = OutputData(filename, prog_average, prog_error, M, l, d);
	
	cout << "Elapsed Times: " << endl;
	cout << "> GenNeedles = " << t1 << " \u03BCs" << endl;
	cout << "> CountNeedles = " << t2 << " \u03BCs" << endl;
	cout << "> BlockAverage = " << t3 << " \u03BCs" << endl;
	cout << "> ConvertToPi = " << t4 << " \u03BCs" << endl;
	cout << "> CopyData = " << t5 << " \u03BCs" << endl;
	cout << "> BlockProgAverage = " << t6 << " \u03BCs" << endl;
	cout << "> OutputData = " << ot1 << " \u03BCs" << endl;

	rnd.SaveSeed();
	return 0;
}
