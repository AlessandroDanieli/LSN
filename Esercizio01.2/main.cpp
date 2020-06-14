#include "random.h"
#include "stat.h"
#include "auxiliary.h"
#include <chrono>

using namespace std;
 
int main (int argc, char *argv[]) {
	Random rnd;
	rnd.RandGenSetup();
	
	/* Parameters */
	unsigned int n_set = 4; // Number of sets represented.
	unsigned int M = 5000000; // Number of measures
	unsigned int L[n_set] = {1, 2, 10, 100}; // Block Lenght
	unsigned int N = 10000;
	string output_file = "dice.out";
	
	double lambda = 1, mean_lorentz = 0, gamma_lorentz = 1;
	vector<double> *dice_std, *dice_exp, *dice_lorentz; // 4 vector<double> each
	dice_std = new vector<double>[n_set];
	dice_exp = new vector<double>[n_set];
	dice_lorentz = new vector<double>[n_set];

	Stat data_unif, data_exp, data_lorentz;
	double t1 = RollDiceContinuous(&rnd, data_unif, data_exp, data_lorentz, M, lambda, mean_lorentz, gamma_lorentz);
	data_unif.SetSampleSize(M); 
	data_exp.SetSampleSize(M); 
	data_lorentz.SetSampleSize(M);
	double t2 = Analyse(data_unif, data_exp, data_lorentz, dice_std, dice_exp, dice_lorentz, n_set, L, N);
	//double ot1 = OutputData("raw_data.out", data_unif, data_exp, data_lorentz, M);
	double ot2 = OutputAnalysis(output_file, dice_std, dice_exp, dice_lorentz, M, N, n_set, L);

	cout << "Elapsed Times: " << endl;
	cout << "> RollDice = " << t1 << " \u03BCs" << endl;
	cout << "> Analyse = " << t2 << " \u03BCs" << endl;
	//cout << "> OutputData = " << ot1 << " \u03BCs" << endl;
	cout << "> OutputAnalysis = " << ot2 << " \u03BCs" << endl;
		
	rnd.SaveSeed();
	return 0;
}
