#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "vector.h"
using namespace std;

Vector :: Vector() {}

Vector :: Vector(unsigned int _d) {
	d = _d;
	V = new double[d];
}

Vector :: Vector(const Vector& x) {
	d = x.GetDim();
	V = new double[d];
	for (unsigned int i = 0; i < d; i++) {
		V[i] = x[i];
	}
	norm = x.norm;
}

Vector :: Vector(const vector<double> x) {
	d = x.size();
	V = new double[d];
	for (unsigned int i = 0; i < d; i++) {
		V[i] = x[i];
	} 
}

/*
// Angular coordinates used as input (identifying the vector with a point in space) || Standard Spherical (N = 3) generalized
Vector :: Vector(double radius, const vector<double> angles) {
	if (radius < 0) {
		cerr << "Error [Vector :: Vector(double, vector<double>)]: negative input radius" << endl;
		exit(-1);
	}
	if (radius == 0) {
		d = angles.size()+1;
		V = new double[d];
		for (unsigned int i = 0; i < d; i++) {
			V[i] = 0;
		}
	}
	else {
		d = angles.size()+1;
		V = new double[d];
		for (unsigned int i = 0; i < d; i++) {
			V[i] = radius;
		}
		
		V[d-1] *= cos(angles[0]); // Last coordinate is simply r*cos(theta_(N-1))
		double s = sin(angles[0]);
		for (unsigned int i = d-2; i >= 1; i--) {
			V[i] *= (s * cos(angles[d-1-i]));
			s *= sin(angles[d-1-i]);
		}
		V[0] *= s;
		
		// *** Swap X1, X2 in order to obtain standard spherical coordinates (generalized to ND): X1 = rsin(theta)cos(phi) for N = 3
		double temp = V[0];
		V[0] = V[1];
		V[1] = temp;
	}
}
*/

Vector :: ~Vector() {
	if (V != nullptr) delete[] V;
}


//## Set Methods ##//

void Vector :: FillValue(unsigned int n, double Value) {
	if (V != nullptr) delete[] V;
	d = n;
	V = new double[d];
	for (unsigned int i = 0; i < d; i++) {
		V[i] = Value;
	}
	norm = -1;
	return;
}

void Vector :: FillValue(double Value) {
	for (unsigned int i = 0; i < d; i++) {
		V[i] = Value;
	}
	norm = -1;
	return;
}

void Vector :: Set(Vector& x) {
	V = x.V;
	d = x.d;
	norm = x.norm;
}

void Vector :: Set(double radius, vector<double> angles) {
	if (radius < 0) {
		cerr << "Error [Vector :: Set(double, vector<double>)]: negative input radius" << endl;
		exit(-1);
	}
	if (angles.size()+1 != d) {
		d = angles.size()+1;
		if (V != nullptr) delete[] V;
		V = new double[d];
	}
	if (radius == 0) {
		for (unsigned int i = 0; i < d; i++) {
			V[i] = 0;
		}
	}
	else {
		for (unsigned int i = 0; i < d; i++) {
			V[i] = radius;
		}
		
		V[d-1] *= cos(angles[0]); // Last coordinate is simply r*cos(theta_(N-1))
		double s = sin(angles[0]);
		for (unsigned int i = d-2; i >= 1; i--) {
			V[i] *= (s * cos(angles[d-1-i]));
			s *= sin(angles[d-1-i]);
		}
		V[0] *= s;
		
		// *** Swap X1, X2 in order to obtain standard spherical coordinates (generalized to ND): X1 = rsin(theta)cos(phi) for N = 3
		double temp = V[0];
		V[0] = V[1];
		V[1] = temp;
	}
	norm = -1;
	return;
}

void Vector :: SetComponent(unsigned int i, double value) {
	if (i < d) V[i] = value;
	else {
		cerr << "Error [Vector :: SetComponent(unsigned int, double)]: index out of Vector bounds" << endl;
		exit(-1);
	}
	norm = -1;
	return;
}

void Vector :: Append(double value) { // Not Efficient!
	if (V != nullptr) {
		double* X = new double[d+1];
		for (unsigned int i = 0; i < d; i++) X[i] = V[i];
		X[d] = value;
		d++;
		delete[] V;
		V = X;
		norm = -1;
	}
	else {
		V = new double[1];
		V[0] = value;
		d = 1;
		norm = value; // Norm == X for a 1D vector
	}
	return;
}

//## Copy Methods ##//

void Vector :: Copy(vector<double> x) {
	if (x.size() == d) {
		for (unsigned int i = 0; i < d; i++) {
			V[i] = x[i];
		}
		norm = -1;
	}
	else {
		if (V != nullptr) delete[] V;
		d = x.size();
		V = new double[d];
		for (unsigned int i = 0; i < d; i++) {
			V[i] = x[i];
		}
		norm = -1;
	}
	return;
}

void Vector :: Copy(Vector x) {
	if (x.d == d) {
		for (unsigned int i = 0; i < d; i++) {
			V[i] = x[i];
		}
	}
	else {
		if (V != nullptr) delete[] V;
		d = x.d;
		V = new double[d];
		for (unsigned int i = 0; i < d; i++) {
			V[i] = x[i];
		}	
	}
	norm = x.norm;
	return;
}


//## Get Methods ##//

double Vector :: GetComponent(unsigned int i) const {
	if (i < d) return V[i];
	else {
		cerr << "Error [Vector :: GetComponent(unsigned int)]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

double Vector :: GetThetaSph() {
	if (d != 3) {
		cerr << "Error [Vector :: GetThetaSph()]: dimension != 3" << endl;
		exit(-1);
	}
	else return acos(V[2]/Norm());
}

double Vector :: GetPhiSph() {
	if (d != 3) {
		cerr << "Error [Vector :: GetPhiSph()]: dimension != 3" << endl;
		exit(-1);
	}
	else return atan2(V[1], V[0]);
}

//## Vector operators ##//

double Vector :: operator[](unsigned int i) const {
	if (i < d) return V[i];
	else {
		cerr << "Error [double Vector :: operator[](unsigned int)]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

double& Vector :: operator[](unsigned int i) {
	if (i < d) return V[i];
	else {
		cerr << "Error [double& Vector :: operator[](unsigned int)]: index out of Vector bounds" << endl;
		exit(-1);
	}
}

bool Vector :: operator==(const Vector& w) const {
	if (d != w.d) return false;
	for (unsigned int i = 0; i < d; i++) {
		if (V[i] != w[i]) return false;
	}
	return true;
}

Vector& Vector :: operator=(const Vector& w) {
	if (V != nullptr) delete[] V;
	d = w.d;
	norm = w.norm;
	V = new double[d];
	for (unsigned int i = 0; i < d; i++) V[i] = w[i];
	return *this; 
}

Vector Vector :: operator+(const Vector& w) const {
	if (d == w.d) {
		Vector sum(w);
		for (unsigned int i = 0; i < d; i++) {
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
	for (unsigned int i = 0; i < d; i++) {
		V[i] = -V[i];
	}
}

Vector Vector :: operator-(const Vector& w) const {
	if (d == w.d) {
		Vector sum(*this);
		for (unsigned int i = 0; i < d; i++) {
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
	if (d == w.d) {
		for (unsigned int i = 0; i < d; i++) {
			scalar_product += V[i] * w[i];
		}
		return scalar_product;
	}
	else {
		cerr << "Error [Vector :: operator*(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
}

Vector Vector :: operator*(double scalar) {
	Vector s(*this);
	for (unsigned int i = 0; i < d; i++) {
		s[i] *= scalar;
	}
	norm *= abs(scalar);
	return s;
}

void Vector :: operator*=(double scalar) {
	*this = *this * scalar;
}

Vector Vector :: operator/(double scalar) {
	if (scalar != 0) {
		Vector s(*this);
		for (unsigned int i = 0; i < d; i++) {
			s[i] /= scalar;
		}
		norm /= abs(scalar);
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

double Vector :: Norm() {
	if (norm == -1) {
		norm = 0;
		for (unsigned int i = 0; i < d; i++) norm += V[i] * V[i];
		norm = sqrt(norm);
		return norm;
	}
	else return norm; // Already calculated
}

double Vector :: Distance(const Vector& Y) const {
	double N = 0;
	for (unsigned int i = 0; i < d; i++) N += (V[i]-Y[i])*(V[i]-Y[i]);
	return sqrt(N);
}

double Vector :: Distance2(const Vector& Y) const { // Distance squared
	double N = 0;
	for (unsigned int i = 0; i < d; i++) N += (V[i]-Y[i])*(V[i]-Y[i]);
	return N;
}

void Vector :: Normalize() {
	double N = Norm();
	if (N == 0) {
		cerr << "Error [Vector :: Normalize]: Norm = 0" << endl;
		exit(-1);
	}
	for (unsigned int i = 0; i < d; i++) V[i] /= N;
	norm = 1;
	return;
}

// Standard n-dim angular coordinates (n-1 angles), reducing to standard spherical for N = 3 [notice that x1 = r*sin(theta)cos(phi)]
// x1 = rsin(t1)...COS(tN)  // x2 = rsin(t1)...SIN(tN)  // x3 = rsin(t1)...COS(t(N-1)) ... // xN = rCOS(t1)
// theta_1, ..., theta_(N-2) -> [0, pi]     while     theta_(N-1) -> [-pi, pi]
vector<double> Vector :: Angles(const string& opt) {
	if (d < 2) {
		cerr << "Error [Vector :: Angles()]: invalid parameters. Dimensions < 2" << endl;
		exit(-1);
	}
	else {
		vector<double> angles;
		double r = Norm();
		
		if (r == 0) { // angles === 0
			for (unsigned int i = 0; i < d-1; i++) {
				angles.push_back(0);
			}
			return angles; 
		}
		else {
			// First Angle
			angles.push_back( acos(V[d-1]/r) ); 
			
			double s = 1, h;
			unsigned int i; 
			for (i = 1; i < d-2; i++) { // The last angle is calculated after this cycle
			
				s *= sin(angles[i-1]); // Progressively multiplying sin(theta_k)
				if (s == 0) { // Remaining angles are arbitrarily set to 0, since they do not contribute
					for (unsigned int k = i; k < d-1; k++) {
						angles.push_back(0);
					}
					break;
				}
				
				else {
					// Manages a floating-point error causing |V[n]/(r*s)| > 1 -> the angle is equal to 0 or PI according to the sign
					if (V[d-1-i] == 0) angles.push_back(M_PI/2.0);
					else {
						h = V[d-1-i]/(r*s);
						if (abs(h) > 1) { // Manages floating point error
							if (h > 0) angles.push_back(0);
							else angles.push_back(M_PI);
						}
						else angles.push_back( acos(h) );
					}
				}
				
			}
			
			// Last angle to be computed: cycle left incomplete (no s == 0 occurred); otherwise already set = 0
			if (i == d-2) {
				angles.push_back(atan2(V[1], V[0])); // ~ atan2(x2/x1) = atan2(y/x)
			}
		}
		
		if (opt == "" || opt == "rad") return angles;
		else if (opt == "deg") {
			for (unsigned int i = 0; i < d-1; i++) {
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
	if (d != w.d) {
		cerr << "Error [Vector :: Parallel(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
	
	if (d != 3) {
		cerr << "Error [Vector :: Parallel(const Vector&)]: dimensions != 3" << endl;
		exit(-1);
	}
	Vector t(*this);
	t.Normalize();
	return t*(w*t);
}

Vector Vector :: Normal(const Vector& w) const {
	if (d != w.d) {
		cerr << "Error [Vector :: Normal(const Vector&)]: different dimensions" << endl;
		exit(-1);
	}
	if (d != 3) {
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

void Vector :: Print(unsigned int prec) const {
	cout.precision(prec);
	if (d == 1) {
		cout << fixed << "x = " << V[0] << endl;
	}
	else if (d == 2) {
		cout << fixed << "(x, y) = (" << V[0] << ", " << V[1] << ")" << endl;
	}
	else if (d == 3) {
		cout << "(x, y, z) = (" << V[0] << ", " << V[1] << ", " << V[2] << ")" << endl;
	}
	else {
		cout << "(";
		for (unsigned int i = 0; i < d-1; i++) {
			cout << "x" << i+1 << ", ";
		}
		cout << "x" << d << ") = (";
		for (unsigned int i = 0; i <= d-2; i++) {
			cout << fixed << V[i] << ", ";
		}
		cout << fixed << V[d-1] << ")" << endl;
	}
	return;
}

void Vector :: PrintSph(const string& opt) {
	cout << setprecision(5);
	
	if (d == 0) cout << "Empty Vector." << endl;
	else if (d == 1) cout << "x = " << V[0] << endl;
	else {
		vector<double> angle = Angles(opt);
		cout << "Spherical coordinates (r, angles): " << endl;
		cout << "Radius = " << Norm() << endl;
		cout << "(Angles 1-" << d - 1 << ") = (";
		for (unsigned int i = 0; i < d-2; i++) {
			cout << angle[i] << ", ";
		}
		cout << angle[d-2] << ")" << endl << endl;
	}
	return;
}

void Vector :: Clear() {
	if (V != nullptr) delete[] V;
	V = nullptr;
	d = 0;
	norm = -1;
}

