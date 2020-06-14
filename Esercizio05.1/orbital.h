#ifndef __Orbital__
#define __Orbital__
#include "vector.h"
using namespace std;

class Orbital {
protected:
	unsigned int Z, n, l;
	int m;
	const double a0 = 5.29177e-11; // Bohr radius in nm
	bool Reduced = true; // Reduced units
	
public:
	Orbital();
	Orbital(unsigned int _Z, unsigned int _n, unsigned int _l, int _m);
	~Orbital();
	
	void SetZ(unsigned int Z); // Set the atomic number
	void SetQN(unsigned int _n, unsigned int _l, int _m); // Set the quantum numbers
	void SetQN(unsigned int Z, unsigned int _n, unsigned int _l, int _m); // Set the quantum numbers
	void SetReduced(bool b) { Reduced = b; } // Set whether to use reduced units (a0) or not
	
	unsigned int GetZ() const { return Z; }
	unsigned int GetN() const { return n; }
	unsigned int GetL() const { return l; }
	int GetM() const { return m; }
	unsigned int GetA0() const { return a0; }
	bool IsReduced() const { return Reduced; }
	
	// Normalized and Not-Normalized functions (due to the radial component). 
	unsigned int Factorial(unsigned int N) const;
	double Binomial(unsigned int N, unsigned int k) const;
	double Laguerre(double x) const; // L[n-l-1, 2l+1]((2*Z/n) * r)
	double R(double r) const; // Radial WaveFunction (r in reduced units) [Normalized]
	double R_nn(double r) const; // Radial WaveFunction (r in reduced units) [Not Normalized]
	double R2(double r) const; // Square of Radial WaveFunction (r in reduced units) [Normalized]
	double R_nn2(double r) const; // Square of Radial WaveFunction (r in reduced units) [Not Normalized]
	double Y_Re(double theta, double phi) const; // Re(Spherical Harmonic) [Normalized]
	double Y_Im(double theta, double phi) const; // Im(Spherical Harmonic) [Normalized]
	// Vector (theta, phi) // Experimental
	//double Y_Re2_nn(Vector X) const; //Re(Spherical Harmonic)^2 [Not Normalized]
	//double Y_Im2_nn(Vector X) const; //Im(Spherical Harmonic)^2 [Not Normalized]
	//
	double Psi_Re(double r, double theta, double phi) const; // Re(WaveFunction) [Normalized] [Requires c++17]
	double Psi_Im(double r, double theta, double phi) const; // Im(WaveFunction) [Normalized] [Requires c++17]
	double Psi_Re_nn(double r, double theta, double phi) const; // Re(WaveFunction) [Not Normalized] [Requires c++17]
	double Psi_Im_nn(double r, double theta, double phi) const; // Im(WaveFunction) [Not Normalized] [Requires c++17]
	double p(double r, double theta) const; // Probability amplitude [Normalized]
	double p_nn(double r, double theta) const; // Probability amplitude [Not Normalized]
	double pv(Vector X);
	double pv_nn(Vector X);
	
	void PrintQN() const;
};

#endif
