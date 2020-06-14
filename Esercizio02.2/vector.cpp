#include "vector.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

Vector :: Vector() {}

Vector :: Vector(const Vector& x) {
	for (unsigned int i = 0; i < x.GetDim(); i++) {
		V.push_back(x[i]);
	}
}

Vector :: Vector(const vector<double> x) {
	for (unsigned int i = 0; i < x.size(); i++) {
		V.push_back(x[i]);
	} 
}

// Angular coordinates used as input (identifying the vector with a point in space)
Vector :: Vector(double radius, const vector<double> angles) {
	if (radius < 0) {
		cerr << "Error [Vector :: Vector(double, vector<double>)]: negative input radius" << endl;
		exit(-1);
	}
	if (radius == 0) {
		for (unsigned int i = 0; i < angles.size()+1; i++) {
			V.push_back(0);
		}
	}
	else {
		for (unsigned int i = 0; i < angles.size()+1; i++) {
			V.push_back(radius);
		}
		
		V[angles.size()] *= cos(angles[0]); // Last coordinate is simply r*cos(theta_(N-1))
		double s = sin(angles[0]);
		for (unsigned int i = angles.size()-1; i >= 1; i--) {
			V[i] *= (s * cos(angles[angles.size()-i]));
			s *= sin(angles[angles.size()-i]);
		}
		V[0] *= s;
		
		// Swap X1, X2 in order to obtain standard spherical coordinates (generalized to ND): X1 = rsin(theta)cos(phi) for N = 3
		double temp = V[0];
		V[0] = V[1];
		V[1] = temp;
	}
}

Vector :: ~Vector() {
	V.clear();
}

//## Set Methods ##//

void Vector :: FillValue(unsigned int n, double Value) {
	V.clear();
	for (unsigned int i = 0; i < n; i++) {
		V.push_back(Value);
	}
	return;
}

void Vector :: FillValue(double Value) {
	for (unsigned int i = 0; i < V.size(); i++) {
		V[i] = Value;
	}
	return;
}

void Vector :: Set(Vector x) {
	V = x.V;
}

void Vector :: Set(vector<double> x) {
	if (x.size() == V.size()) {
		for (unsigned int i = 0; i < x.size(); i++) {
			V[i] = x[i];
		}
	}
	else {
		V.clear();
		for (unsigned int i = 0; i < x.size(); i++) {
			V.push_back(x[i]);
		}	
	}
	return;
}

void Vector :: Set(double radius, vector<double> angles) {
	if (radius < 0) {
		cerr << "Error [Vector :: Set(double, vector<double>)]: negative input radius" << endl;
		exit(-1);
	}
	V.clear();
	if (radius == 0) {
		for (unsigned int i = 0; i < angles.size()+1; i++) {
			V.push_back(0);
		}
	}
	else {
		for (unsigned int i = 0; i < angles.size()+1; i++) {
			V.push_back(radius);
		}
		
		V[angles.size()] *= cos(angles[0]); // Last coordinate is simply r*cos(theta_(N-1))
		double s = sin(angles[0]);
		for (unsigned int i = angles.size()-1; i >= 1; i--) {
			V[i] *= (s * cos(angles[angles.size()-i]));
			s *= sin(angles[angles.size()-i]);
		}
		V[0] *= s;
		// Swap X1, X2 in order to obtain standard spherical coordinates (generalized to ND): X1 = rsin(theta)cos(phi) for N = 3
		double temp = V[0];
		V[0] = V[1];
		V[1] = temp;
	}
	return;
}

void Vector :: SetComponent(unsigned int i, double x) {
	if (i >= 0 && i < V.size()) V[i] = x;
	else {
		cerr << "Error [Vector :: SetComponent]: index out of Vector bounds" << endl;
		exit(-1);
	}
	return;
}

void Vector :: Append(double value) {
	V.push_back(value);
	return;
}

//## Copy Methods ##//

void Vector :: Copy(vector<double>& P) {
	if (P.size() == 0) {
		cerr << "Error [Vector :: Copy(vector<double>)]: empty input coordinate set" << endl;
		exit(-1);
	}
	else {
		V.clear();
		for (unsigned int i = 0; i < P.size(); i++) {
			V.push_back(P[i]);
		}
	}
	return;
}

void Vector :: Copy(Vector x) {
	V.clear();
	for (unsigned int i = 0; i < x.GetDim(); i++) {
		V.push_back(x[i]);
	}
	return;
}

//## Get Methods ##//

double Vector :: GetComponent(unsigned int i) const {
	if (i >= 0 && i < V.size()) return V[i];
	else {
		cerr << "Error [Vector :: GetComponent]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

//## Vector operators ##//

double Vector :: operator[](unsigned int i) const {
	if (i >= 0 && i < V.size()) return V[i];
	else {
		cerr << "Error [double Vector :: operator[](unsigned int)]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

double& Vector :: operator[](unsigned int i) {
	if (i >= 0 && i < V.size()) return V[i];
	else {
		cerr << "Error [double& Vector :: operator[](unsigned int)]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

bool Vector :: operator==(const Vector& w) const {
	if (V.size() != w.GetDim()) return false;
	for (unsigned int i = 0; i < V.size(); i++) {
		if (V[i] != w[i]) return false;
	}
	return true;
}

Vector& Vector :: operator=(const Vector& w) {
	V.clear();
	for (unsigned int i = 0; i < w.GetDim(); i++) V.push_back(w[i]);
	return *this; 
}

Vector Vector :: operator+(const Vector& w) const {
	if (V.size() == w.GetDim()) {
		Vector sum(w);
		for (unsigned int i = 0; i < V.size(); i++) {
			sum[i] += V[i];
		}
		return sum;
	}
	
	else {
		cerr << "Error [Vector :: operator+(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
}

void Vector :: operator+=(const Vector& w) {
	*this = *this + w;
}

void Vector :: operator-()  {
	for (unsigned int i = 0; i < V.size(); i++) {
		V[i] = -V[i];
	}
}

Vector Vector :: operator-(const Vector& w) const {
	if (V.size() == w.GetDim()) {
		Vector sum(*this);
		for (unsigned int i = 0; i < V.size(); i++) {
			sum[i] -= w[i];
		}
		return sum;
	}
	else {
		cerr << "Error [Vector :: operator-(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
}

void Vector :: operator-=(const Vector& w) {
	*this = *this - w;
}

//## Scalar product and scalar operators ##//

double Vector :: operator*(const Vector& w) const { //Scalar Product
	double scalar_product = 0;
	if (V.size() == w.GetDim()) {
		for (unsigned int i = 0; i < V.size(); i++) {
			scalar_product += V[i] * w[i];
		}
		return scalar_product;
	}
	else {
		cerr << "Error [Vector :: operator*(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
}

Vector Vector :: operator*(double scalar) const {
	Vector s(*this);
	for (unsigned int i = 0; i < V.size(); i++) {
		s[i] *= scalar;
	}
	return s;
}

void Vector :: operator*=(double scalar) {
	*this = *this * scalar;
}

Vector Vector :: operator/(double scalar) const {
	if (scalar != 0) {
		Vector s(*this);
		for (unsigned int i = 0; i < V.size(); i++) {
			s[i] /= scalar;
		}
		return s;
	}
	else {
		cerr << "Error [Vector :: operator/(double)]: division by zero" << endl;
		exit(-1);
	}
}

void Vector :: operator/=(double scalar) {
	*this = *this / scalar;
}

double Vector :: Norm() const {
	double N = 0;
	for (unsigned int i = 0; i < V.size(); i++) N += V[i] * V[i];
	return sqrt(N);
}

double Vector :: Norm(const Vector& Y) const {
	double N = 0;
	for (unsigned int i = 0; i < V.size(); i++) N += (V[i]-Y[i])*(V[i]-Y[i]);
	return sqrt(N);
}

void Vector :: Normalize() {
	double N = Norm();
	if (N == 0) {
		cerr << "Error [Vector :: Normalize]: Norm = 0" << endl;
		exit(-1);
	}
	for (unsigned int i = 0; i < V.size(); i++) V[i] /= N;
	return;
}

// Standard n-dim angular coordinates (n-1 angles), reducing to standard spherical for N = 3 (note x1 = rsin(theta)cos(phi)
// x1 = rsin(t1)...COS(tN)  // x2 = rsin(t1)...SIN(tN)  // x3 = rsin(t1)...COS(t(N-1)) ... // xN = rCOS(t1)
// theta_1, ... theta_(N-2) -> [0, pi] while theta_(N-1) -> [-pi, pi]
vector<double> Vector :: Angles(const string& opt) const {
	if (V.size() < 2) {
		cerr << "Error [Vector :: Angles()]: invalid parameters. Dimensions < 2" << endl;
		exit(-1);
	}
	else {
		vector<double> angles;
		double r = Norm();
		
		if (r == 0) {
			for (unsigned int i = 0; i < V.size()-1; i++) {
				angles.push_back(0);
			}
			return angles; // angles === 0
		}
		else {
			angles.push_back( acos(V[V.size()-1]/r) );
			
			double s = 1;
			unsigned int i;
			for (i = 1; i < V.size()-2; i++) {
				s *= sin(angles[i-1]); // Progressively multiplying sin(theta_k)
				if (s == 0) { // Remaining angles are arbitrarily set to 0, since they do not contribute
					for (unsigned int k = i; k < V.size()-2; k++) {
						angles[k] = 0;
					}
					break;
				}
				else {
					angles[i] = acos(V[V.size()-1-i]/(r*s));
				}
			}
			// Last angle to be computed: cycle left incomplete (no s == 0 occurred); otherwise already set = 0
			if (i == V.size()-2) { 
				angles[V.size()-2] = atan2(V[1], V[0]);
			}
		}
		
		if (opt == "" || opt == "rad") return angles;
		else if (opt == "deg") {
			for (unsigned int i = 0; i < V.size()-1; i++) {
				angles[i] *= 180./M_PI;
			}
			return angles;
		}
		else {
			cerr << "Error [Vector :: Angles()]: invalid option (possible: rad (default), deg)" << endl;
			exit(-1);
		}	
	}
}

Vector Vector :: Parallel(const Vector& w) const {
	if (V.size() != w.GetDim()) {
		cerr << "Error [Vector :: Parallel(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
	
	if (V.size() != 3) {
		cerr << "Error [Vector :: Parallel(const Vector&)]: dimensions != 3" << endl;
		exit(-1);
	}
	Vector t(*this);
	t.Normalize();
	return t*(w*t);
}

Vector Vector :: Normal(const Vector& w) const {
	if (V.size() != w.GetDim()) {
		cerr << "Error [Vector :: Normal(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
	if (V.size() != 3) {
		cerr << "Error: [Vector :: Normal(const Vector&)] != 3" << endl;
		exit(-1);
	}
	return w - this->Parallel(w);
}


void Vector :: Swap(unsigned int p, unsigned int q) {
	double temp = V[p];
	V[p] = V[q];
	V[q] = temp;
}

void Vector :: Print() const {
	if (V.size() == 1) {
		cout << "x = " << V[0] << endl;
	}
	else if (V.size() == 2) {
		cout << "(x, y) = (" << V[0] << ", " << V[1] << ")" << endl;
	}
	else if (V.size() == 3) {
		cout << "(x, y, z) = (" << V[0] << ", " << V[1] << ", " << V[2] << ")" << endl;
	}
	else {
		cout << "(";
		for (unsigned int i = 0; i <= V.size()-2; i++) {
			cout << V[i] << ", ";
		}
		cout << V[V.size()-1] << ")" << endl << endl;
	}
	return;
}

void Vector :: PrintSph(const string& opt) const {
	cout << setprecision(5);
	
	if (V.size() == 0) cout << "Empty Vector." << endl;
	else if (V.size() == 1) cout << "x = " << V[0] << endl;
	else {
		vector<double> angle = Angles(opt);
		cout << "Spherical coordinates (r, angles): " << endl;
		cout << "Radius = " << Norm() << endl;
		cout << "(Angles 1-" << V.size() - 1 << ") = (";
		for (unsigned int i = 0; i < V.size()-2; i++) {
			cout << angle[i] << ", ";
		}
		cout << angle[V.size()-2] << ")" << endl << endl;
	}
	return;
}
