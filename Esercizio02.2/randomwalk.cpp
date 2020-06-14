#include <iostream>
#include <cmath>
#include "randomwalk.h"

using namespace std;

Randomwalk :: Randomwalk(): Rnd(nullptr) {}

Randomwalk :: Randomwalk(Random* r): Rnd(r) {}

Randomwalk :: Randomwalk(Random* r, vector<double> x): Rnd(r) {
	if (x.size() == 0) {
		cerr << "Error [Randomwalk :: RandomMove(Random*, vector<double>)]: empty input coordinate set" << endl;
		exit(-1);
	}
	P.Copy(x);
}

Randomwalk :: Randomwalk(Random* r, Vector p): Rnd(r), P(p) {}

Randomwalk :: ~Randomwalk() {
	if (Rnd != nullptr) Rnd->SaveSeed();
	Rnd = nullptr;
}

// Using Random :: Uniform.
void Randomwalk :: RandomMove(double distance, unsigned int k) {
	if (k == 0) return;
	if (distance < 0) {
		cerr << "Error [Randomwalk :: RandomMove(double, unsigned int)]: negative input distance" << endl;
		exit(-1);
	}
	if (P.GetDim() > 1) {
		vector<double> angles; // Angular coordinates
		
		for (unsigned int i = 0; i < P.GetDim()-2; i++) {
			angles.push_back(Rnd->Uniform(0, M_PI)); // from 0 to PI
		}
		angles.push_back(Rnd->Uniform(0, 2.*M_PI)); // from 0 to 2PI
		
		Vector delta(distance, angles), t;

		for (unsigned int j = 0; j < k-1; j++) {
			//angles.clear();
			for (unsigned int i = 0; i < P.GetDim()-2; i++) {
				angles[i] = Rnd->Uniform(0, M_PI); // from 0 to PI
			}
			angles[P.GetDim()-2] = Rnd->Uniform(0, 2.*M_PI); // from 0 to 2PI
			t.Set(distance, angles);
			delta += t;
		}
		P += delta;
		return;
	}
	if (P.GetDim() == 1)  {
		for (unsigned int j = 0; j < k; j++) {
			P[0] += distance * (2*Rnd->Bool() - 1);
		}
		return;
	}
}

// Using Random :: Angle to avoid correlation with the value of PI. Less efficient than the alternative version showed above
/*
void Randomwalk :: RandomMove(double distance) {
	if (distance < 0) {
		cerr << "Error [Randomwalk :: RandomMove(double)]: negative input distance" << endl;
		exit(-1);
	}
	vector<double> angles; // Angular coordinates
	for (unsigned int i = 0; i < P.GetDim()-2; i++) {
		angles.push_back(Rnd->Angle()/2.); // from 0 to PI
	}
	angles.push_back(Rnd->Angle()); // from 0 to 2PI
	
	Vector delta(distance, angles);
	pos += delta;
	return;
}
*/

double Randomwalk :: Distance() const {
	return P.Norm();
}

double Randomwalk :: Distance(const vector<double>& x) {
	if (x.size() != P.GetDim()) {
		cerr << "Error: [Randomwalk :: Distance(const vector<double>&)]: input coordinate set and Position member have different\
		 dimensions" << endl;
		exit(-1);
	}
	Vector X(x);
	return P.Norm(X);
}

double Randomwalk :: Distance(const Vector& x) {
	if (x.GetDim() != P.GetDim()) {
		cerr << "Error: [Randomwalk :: Distance(const Vector&)]: input coordinate set and Position member have different\
		 dimensions" << endl;
		exit(-1);
	}
	return P.Norm(x);
}

double Randomwalk :: Distance(const Randomwalk& x) {
	if (x.GetPosition().GetDim() != P.GetDim()) {
		cerr << "Error: [Randomwalk :: Distance(const Randomwalk&)]: input Randomwalk and (*this) have different Position dimensions"\
		 << endl;
		exit(-1);
	}
	Vector X(x.GetPosition());
	return P.Norm(X);
}

void Randomwalk :: PrintPosition() const {
	cout << "Actual position: " << endl;
	P.Print();
	return;
}


// ### RandomwalkLattice Class ### // // ### RandomwalkLattice Class ### // // ### RandomwalkLattice Class ### //
RandomwalkLattice :: RandomwalkLattice() {
	Rnd = nullptr;
}

RandomwalkLattice :: RandomwalkLattice(Random* r, Vector p, vector<Vector> g) {
	Rnd = r;
	P = p;
	if (P.GetDim() < g.size()) {
		cerr << "Error [RandomwalkLattice :: RandomwalkLattice(Random*, Vector, vector<Vector>]: the number of generators exceeds the space dimension. " << endl;
		exit(-1);
	}
	Gen = g;
}

RandomwalkLattice :: ~RandomwalkLattice() {
	if (Rnd != nullptr) Rnd->SaveSeed();
	Rnd = nullptr;
}

void RandomwalkLattice :: SetGenerators(const vector<Vector>& g) {
	if (P.GetDim() < g.size()) {
		cerr << "Error [RandomwalkLattice :: SetGenerators(vector<Vector>]: the number of generators exceeds the space dimension. " << endl;
		exit(-1);
	}
	Gen = g;
}

void RandomwalkLattice :: RandomMove(const unsigned int k) {
	// Direction, Orientation
	for (unsigned int i = 0; i < k; i++) {
		if (Rnd->Bool()) {
			P += Gen[(unsigned int) Rnd->Uniform(0, Gen.size())];
		}
		else {
			P -= Gen[(unsigned int) Rnd->Uniform(0, Gen.size())];
		}
	}
	return;
}


