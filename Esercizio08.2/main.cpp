#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"
#include "stat.h"
#include "variational.h"
using namespace std;

int main() {
	Random Rnd; 
	Rnd.RandGenSetup();
	
	Variational V(&Rnd);
	V.Input("input.dat");
	V.Minimize();
	V.Calibrate();
	V.Simulate();
	V.OutputSample();
	
	Rnd.SaveSeed();
	return 0;
}



