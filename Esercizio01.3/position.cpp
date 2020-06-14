#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "position.h"

// Angles(): Fix each theta being in [0, 2pi] //
/// /// // // // / / / / / // / / / // / // / //
using namespace std;
Position :: Position() {}

Position :: Position(double x, double y) {
	X.push_back(x);
	X.push_back(y);
}

//L'inserimento avviene sempre utilizzando coordinate cartesiane, ma tali coordinate vengono convertite in sferiche internamente
Position :: Position(double x, double y, double z): Position(x, y) {
	X.push_back(z);
}

Position :: Position(vector<double> P) {
	if (P.size() == 0) {
		cerr << "Error [Position :: Position(vector<double> )]: empty set of coordinates" << endl;
	}
	else {
		for (unsigned int i = 0; i < P.size(); i++) {
			X.push_back(P[i]);
		}
	}
}

Position :: Position(const Position& P) {
	X.clear();
	for (unsigned int i = 0; i < P.GetDim(); i++) {
		X.push_back(P.X[i]);
	}
}

Position :: ~Position() {
	X.clear();
}

Position& Position :: operator=(const Position& P) {
	X.clear();
	for (unsigned int i = 0; i < P.GetDim(); i++) {
		X.push_back(P.X[i]);
	}
	return *this;
}

//## Set Methods ##//
void Position :: Set(double x, double y) {
	if (X.size() != 2) {
		cerr << "Error [Position :: Set(double, double)]: dimensions != 2" << endl;
		exit(-1);
	}
	else {
		X[0] = x;
		X[1] = y;
	}
}
void Position :: Set(double x, double y, double z) {
	if (X.size() != 3) {
		cerr << "Error [Position :: Set(double, double, double)]: dimensions != 3" << endl;
		exit(-1);
	}
	else {
		X[0] = x;
		X[1] = y;
		X[2] = z;
	}
}

void Position :: Set(vector<double> P) {
	if (P.size() != X.size()) {
		cerr << "Error [Position :: Set(vector<double>)]: invalid dimension of input vector" << endl;
		exit(-1);
	}
	else {
		for (unsigned int i = 0; i < X.size(); i++) {
			X[i] = P[i];
		}
	}
	return;
}

//## Get Methods ##/
double Position :: GetR() const { 
	double r = 0;
	for (unsigned int i = 0; i < X.size(); i++) {
		r += X[i]*X[i];
	}
	return sqrt(r); 
}

double Position::GetX() const { return X[0]; }
double Position::GetY() const { 
	if (X.size() < 2) {
		cerr << "Error [Position :: GetY()]: invalid coordinate (dimensions < 2)" << endl;
		exit(-1);
	}
	else {
		return X[1]; 
	}
}
double Position::GetZ() const {
	if (X.size() < 3) {
		cerr << "Error [Position :: GetZ()]: invalid coordinate (dimensions < 3)" << endl;
		exit(-1);
	}
	else {
		return X[2]; 
	}
}
double Position :: GetRho() const { 
	if (X.size() != 3) {
		cerr << "Error [Position :: GetPhi()]: invalid coordinate (dimensions != 3)" << endl;
		exit(-1);
	}
	else return sqrt(X[0]*X[0] + X[1]*X[1]); 
}
double Position :: GetPhi(const string& opt) const { 
	if (X.size() != 3) {
		cerr << "Error [Position :: GetPhi()]: invalid coordinate (dimensions != 3)" << endl;
		exit(-1);
	}
	else if (opt == "" || opt == "rad") return atan2(X[1], X[0]);
	else if (opt == "deg") return atan2(X[1], X[0])*180./M_PI;
	else {
		cerr << "Error [Position :: GetPhi()]: invalid option (possible: rad (default), deg)" << endl;
		exit(-1);
	}
}
double Position :: GetTheta(const string& opt) const {
	if (X.size() != 2 && X.size() != 3) {
		cerr << "Error [Position :: GetTheta()]: invalid coordinate (dimensions != 2 or 3)" << endl;
		exit(-1);
	}
	else if (X.size() == 2) { // Polar 2D
		if (opt == "" || opt == "rad") return GetR() == 0 ? 0 : atan2(X[1], X[0]);
		else if (opt == "deg") return GetR() == 0 ? 0 : atan2(X[1], X[0])*180./M_PI;
		else {
			cerr << "Error [Position :: GetTheta()]: invalid option (possible: rad (default), deg)" << endl;
			exit(-1);
		}
	}
	else { // Spherical 3D
		if (opt == "" || opt == "rad") return GetR() == 0 ? 0 : acos(X[2]/GetR());
		else if (opt == "deg") return GetR() == 0 ? 0 : acos(X[2]/GetR())*180./M_PI;
		else {
			cerr << "Error [Position :: GetTheta()]: invalid option (possible: rad (default), deg)" << endl;
			exit(-1);
		}
	}
}

// Standard n-dim angular coordinates (n-1 angles)
// x1 = rsin(t1)...sin(tN)  // x2 = rsin(t1)...cos(tN)  // x3 = rsin(t1)...cos(t(N-1)) ... // xN = rcos(t1)
// theta_1, ... theta_(N-2) -> [0, pi] while theta_(N-1) -> [-pi, pi]
vector<double> Position :: Angles(const string& opt) const {
	if (X.size() < 2) {
		cerr << "Error [Position :: Angles()]: invalid call. Dimensions < 2" << endl;
		exit(-1);
	}
	else {
		vector<double> angles;
		double r = GetR();
		
		if (r == 0) {
			for (unsigned int i = 0; i < X.size()-1; i++) {
				angles[i] = 0;
			}
			return angles; // angles === 0
		}
		else {
			angles.push_back( acos(X[X.size()-1]/r) );
			
			double s = 1;
			unsigned int i;
			for (i = 1; i < X.size()-2; i++) {
				s *= sin(angles[i-1]); // Progressively multiplying sin(theta_k)
				if (s == 0) { // Remaining angles are arbitrarily set to 0, since they do not contribute
					for (unsigned int k = i; k < X.size()-2; k++) {
						angles[k] = 0;
					}
					break;
				}
				else {
					angles[i] = acos(X[X.size()-1-i]/(r*s));
				}
			}
			if (i == X.size()-2) { // Last angle to be computed: cycle left incomplete (no s == 0 occurred); otherwise already set = 0
				angles[X.size()-2] = atan2(X[0], X[1]);
			}
		}
		
		if (opt == "" || opt == "rad") return angles;
		else if (opt == "deg") {
			for (unsigned int i = 0; i < X.size()-1; i++) {
				angles[i] *= 180./M_PI;
			}
			return angles;
		}
		else {
			cerr << "Error [Position :: Angles()]: invalid option (possible: rad (default), deg)" << endl;
			exit(-1);
		}	
	}
}

//Distanza deall'origine
double Position :: Distance() const {
	double r = 0;
	for (unsigned int i = 0; i < X.size(); i++) {
		r += X[i]*X[i];
	}
	return sqrt(r);
}

//Distanza da una posizione generica P
double Position :: Distance(const Position& P) const {
	if (P.GetDim() != X.size()) {
		cerr << "Error [Position :: Distance(const Position&)]: invalid dimension of input vector" << endl;
		exit(-1);
	}
	double r = 0;
	for (unsigned int i = 0; i < X.size(); i++) {
		r += (X[i]-P.X[i])*(X[i]-P.X[i]);
	}
	return sqrt(r);
}

void Position :: Print() const {
	cout << setprecision(5);
	if (X.size() == 1) {
		cout << "x = " << X[0] << endl;
	}
	else if (X.size() == 2) {
		cout << "(x, y) = (" << GetX() << ", " << GetY() << ")" << endl;
	}
	else if (X.size() == 3) {
		cout << "(x, y, z) = (" << GetX() << ", " << GetY() << ", " << GetZ() << ")" << endl;
	}
	else {
		cout << "(";
		for (unsigned int i = 0; i < X.size()-1; i++) {
			cout << X[i] << ", ";
		}
		cout << X[X.size()-1] << ")" << endl;
	}	
	return;
}

void Position :: PrintPolar(const string& opt) const {
	cout << setprecision(5);
	if (X.size() != 2) {
		cerr << "Error [Position :: PrintPolar()]: dimensions != 2" << endl;
		exit(-1);
	}
	else cout << "(r, theta) = (" << GetR() << ", " << GetTheta(opt) << ")" << endl;
	return;
}

void Position :: PrintCyl(const string& opt) const {
	cout << setprecision(5);
	if (X.size() != 3) {
		cerr << "Error [Position :: PrintCyl()]: dimensions != 3" << endl;
		exit(-1);
	}
	else cout << "(rho, phi, z) = (" << GetRho() << ", " << GetPhi(opt) << ", " << GetZ() << ")" << endl;
	return;
}

void Position::PrintSph(const string& opt) const {
	cout << setprecision(5);
	
	if (X.size() == 3) {
		cout << "(r, phi, theta) = (" << GetR() << ", " << GetPhi(opt) << ", " << GetTheta(opt) << ")" << endl;
	}
	else {
		vector<double> angle = Angles(opt);
		cout << "Radius = " << GetR() << endl;
		cout << "(Angles 1-" << X.size() - 1 << ") = (";
		for (unsigned int i = 0; i < X.size()-2; i++) {
			cout << angle[i] << ", ";
		}
		cout << angle[X.size()-2] << ")" << endl;
	}
	return;
}
