#include <iostream>
#include "function.h"

using namespace std;

Gauss :: Gauss(double m, double s) {
	if (s == 0) {
		cerr << "Error [Gauss :: Gauss(double, double)]: sigma = 0" << endl;
		exit(-1);
	}
	mean = m;
	sigma = s;
	k = sigma*sqrt(2.*M_PI); // Normalized;
}

void Gauss :: Normalize() {
	if (sigma == 0) {
		cerr << "Error [Gauss :: Normalize()]: sigma = 0" << endl;
		exit(-1);
	}
	k = 1./(sigma*sqrt(2.*M_PI));
	return;
}

double Gauss :: Eval(double x) const {
	if (sigma == 0) {
		cerr << "Error [Gauss :: Eval(double)]: sigma = 0" << endl;
		exit(-1);
	}
 	return k*exp(-(x-mean)*(x-mean)/(2.*sigma*sigma));
}

Lorentz :: Lorentz(double m, double g) {
	if (g == 0) {
		cerr << "Error [Lorentz :: Lorentz(double, double)]: gamma = 0" << endl;
		exit(-1);
	}
	mean = m;
	gamma = g;
	k = gamma*sqrt(2.*M_PI); // Normalized;
}

void Lorentz :: Normalize() {
	if (gamma == 0) {
		cerr << "Error [Lorentz :: Normalize()]: gamma = 0" << endl;
		exit(-1);
	}
	k = gamma/M_PI;
}

