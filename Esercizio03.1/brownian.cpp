#include "brownian.h"
#include "random.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

Brownian :: Brownian(): rnd(nullptr), w(0), t(0), mean(0), sigma(1) {}

Brownian :: Brownian(Random* r): w(0), t(0), mean(0), sigma(1) {
	if (r == nullptr) {
		cerr << "Error [Brownian :: Brownian(Random*)]: Random* argument is a nullptr" << endl;
		exit(-1);
	}
	rnd = r;
}

Brownian :: Brownian(Random* r, double w0, double t0, double _mean, double _sigma): rnd(r), w(w0), t(t0), mean(_mean), sigma(_sigma) {}

Brownian :: ~Brownian() {
	if (rnd != nullptr) rnd->SaveSeed();
}

void Brownian :: SetRandGen(Random* r) {
	if (r == nullptr) {
		cerr << "Error [Brownian :: SetRandGen(Random*)]: Random* argument is a nullptr" << endl;
		exit(-1);
	}
	rnd = r;
}

void Brownian :: Move() { // Assuming Delta(t) = 1
	t++;
	w += rnd->Gauss(mean, sigma);
	return;
}
void Brownian :: Move(double delta_t) { // variance = Delta(t)
	t += delta_t;
	w += mean*delta_t + sigma*sqrt(delta_t)*rnd->Gauss();
	return;
}



GBM :: GBM(): Brownian() {}
GBM :: GBM(Random* r): Brownian(r) {}
GBM :: GBM(Random* r, double mean, double sigma, double S0): Brownian(r, 0, 0, mean, sigma), S(S0) {}
GBM :: ~GBM() {}

void GBM :: Evolve(double delta_t, unsigned int steps) {
	if (delta_t < 0) {
		cerr << "Error [GBM :: Evolve(double, unsigned int)]: delta_t must be positive" << endl;
		exit(-1);
	}
	for (unsigned int n = 0; n < steps; n++) {
		t += delta_t;
		S *= exp( (mean-sigma*sigma/2.)*delta_t + sigma*sqrt(delta_t)*rnd->Gauss() );
	}
	return;
}
