#include <cmath>
#include <iostream>
#include "orbital.h"
#include "vector.h"
using namespace std;

Orbital :: Orbital(): Z(1), n(1), l(0), m(0), Reduced(true) {} // Hydrogen atom by default

Orbital :: Orbital(unsigned int _Z, unsigned int _n, unsigned int _l, int _m) {
	if (_Z == 0) {
		cerr << "Error [Orbital :: Orbital(unsigned int, unsigned int, unsigned int)]: Invalid input Z" << endl;
		exit(-1);
	}
	if (_n == 0 || _l >= _n || abs(_m) > _l) {
		cerr << "Error [Orbital :: Orbital(unsigned int, unsigned int, unsigned int, int)]: Invalid quantum numbers" << endl;
		exit(-1);
	}
	else {
		Z = _Z;
		n = _n;
		l = _l;
		m = _m;
		Reduced = true;
	}
}

Orbital :: ~Orbital() {}

void Orbital :: SetZ(unsigned int _Z) {
	if (_Z == 0) {
		cerr << "Error [Orbital :: SetZ(unsigned int)]: Invalid input Z" << endl;
		exit(-1);
	}
	else Z = _Z;
	return;
}

void Orbital :: SetQN(unsigned int _n, unsigned int _l, int _m) {
	if (_n == 0 || _l >= _n || abs(_m) > _l) {
		cerr << "Error [Orbital :: SetQN(unsigned int, unsigned int, int)]: Invalid quantum numbers" << endl;
		exit(-1);
	}
	else {
		n = _n;
		l = _l;
		m = _m;
	}
	return;
}

void Orbital :: SetQN(unsigned int _Z, unsigned int _n, unsigned int _l, int _m) {
	SetZ(_Z);
	SetQN(_n, _l, _m);
	return;
}

unsigned int Orbital :: Factorial(unsigned int N) const {
	if (N == 0 || N == 1) return 1;
	else {
		for (unsigned int i = N-1; i >= 2; i--) N *= i;
		return N;
	}
}

double Orbital :: Binomial(unsigned int N, unsigned int k) const {
	if (N < k) {
		cerr << "Error [Orbital :: Binomial]: N > k" << endl;
		exit(-1);
	}
	else if (N == k || (N == 1 && k == 0)) return 1;
	else {
		for (unsigned int i = N-1; i >= N-k+1; i--) N *= i;
		return N/(double)Factorial(k);
	}
}

double Orbital :: Laguerre(double x) const {
	unsigned int N = n-l-1, a = 2*l+1;
	double L = 0;
	for (unsigned int k = 0; k <= N; k++) {
		if (k % 2 == 0) L += Binomial(N+a, k+a)/((double)Factorial(k))*pow(x, k);
		else L -= Binomial(N+a, k+a)/((double)Factorial(k))*pow(x, k);
	}
	return L;
}

double Orbital :: R(double r) const {
	if (!Reduced) r /= a0;
	double rho = (2*Z/(double)n) * r;
	return pow(2*Z/(double)n, 3./2.)*sqrt(Factorial(n-l-1)/( (double) 2*n*Factorial(n+l) ))*pow(rho, l)*exp(-rho/2.)*Laguerre(rho);
}

double Orbital :: R2(double r) const { // Just R(r) * R(r)
	if (!Reduced) r /= a0;
	double rho = (2*Z/(double)n) * r;
	double L = Laguerre(rho);
	return pow(2*Z/(double)n, 3.)*Factorial(n-l-1) / ( (double)2*n*Factorial(n+l) ) *pow(rho, 2*l)*exp(-rho)*L*L;
}

double Orbital :: R_nn(double r) const {
	if (!Reduced) r /= a0;
	double rho = (2*Z/(double)n) * r;
	return pow(rho, l)*exp(-rho/2.)*Laguerre(rho);
}

double Orbital :: R_nn2(double r) const { // Just R_nn(r) * R_nn(r)
	if (!Reduced) r /= a0;
	double rho = (2*Z/(double)n) * r;
	double L = Laguerre(rho);
	return pow(rho, 2*l)*exp(-rho)*L*L;
}

double Orbital :: Y_Re(double theta, double phi) const {
	if (m >= 0) return sph_legendre(l, m, theta)*cos(m*phi); // Associated Legendre Polynomial, m >= 0
	else return (1-2*(-m%2))*Factorial(l+m)/(double)Factorial(l-m)*sph_legendre(l, -m, theta)*cos(m*phi);
}
double Orbital :: Y_Im(double theta, double phi) const {
	if (m >= 0) return sph_legendre(l, m, theta)*sin(m*phi); // Associated Legendre Polynomial, m >= 0
	else return (1-2*(-m%2))*Factorial(l+m)/(double)Factorial(l-m)*sph_legendre(l, -m, theta)*sin(m*phi);
}

/*
double Orbital :: Y_Re2_nn(Vector X) const {
	double SH;
	if (m >= 0) SH = sph_legendre(l, m, X[0])*cos(m*X[1]); // Associated Legendre Polynomial, m >= 0
	else SH = sph_legendre(l, -m, X[0])*cos(m*X[1]);
	return SH*SH;
}
double Orbital :: Y_Im2_nn(Vector X) const {
	double SH;
	if (m >= 0) SH = sph_legendre(l, m, X[0])*sin(m*X[1]); // Associated Legendre Polynomial, m >= 0
	else SH = sph_legendre(l, -m, X[0])*sin(m*X[1]);
	return SH*SH;
}
*/
double Orbital :: Psi_Re(double r, double theta, double phi) const {
	return R(r)*Y_Re(theta, phi);
}

double Orbital :: Psi_Im(double r, double theta, double phi) const {
	return R(r)*Y_Im(theta, phi);
}

double Orbital :: Psi_Re_nn(double r, double theta, double phi) const {
	return R_nn(r)*Y_Re(theta, phi);
}

double Orbital :: Psi_Im_nn(double r, double theta, double phi) const {
	return R_nn(r)*Y_Im(theta, phi);
}

double Orbital :: p(double r, double theta) const {
	double SH;
	if (m >= 0) SH = sph_legendre(l, m, theta); // Associated Legendre Polynomial, m >= 0
	else SH = Factorial(l+m)/(double)Factorial(l-m)*sph_legendre(l, -m, theta);
	return R2(r)*SH*SH;
}

double Orbital :: p_nn(double r, double theta) const {
	double SH;
	if (m >= 0) SH = sph_legendre(l, m, theta); // Associated Legendre Polynomial, m >= 0
	else SH = sph_legendre(l, -m, theta); // Normalization of P(l, -m) discarded
	return R_nn2(r)*SH*SH;
}

double Orbital :: pv(Vector X) {
	if (X.GetDim() != 3) {
		cerr << "Error [Orbital :: pv(Vector)]: the input Vector has dimension != 3" << endl;
		exit(-1);
	}
	else {
		double SH;
		if (m >= 0) SH = sph_legendre(l, m, X.GetThetaSph()); // Associated Legendre Polynomial, m >= 0
		else SH = Factorial(l+m)/(double)Factorial(l-m)*sph_legendre(l, -m, X.GetThetaSph());
		return R2(X.Norm())*SH*SH;
	}
}

double Orbital :: pv_nn(Vector X) { // Useful when using the Metropolis Algorithm
	if (X.GetDim() != 3) {
		cerr << "Error [Orbital :: pv_nn(Vector)]: the input Vector has dimension != 3" << endl;
		exit(-1);
	}
	else {
		double SH;
		if (m >= 0) SH = sph_legendre(l, m, X.GetThetaSph()); // Associated Legendre Polynomial, m >= 0
		else SH = sph_legendre(l, -m, X.GetThetaSph()); // Normalization of P(l, -m) discarded
		return R_nn2(X.Norm())*SH*SH;
	}
}

void Orbital :: PrintQN() const {
	cout << "Orbital Quantum Numbers: (Z, n, l, ml) = (" << Z << ", " << n << ", " << l << ", " << m << ")" << endl;
	return;
}
