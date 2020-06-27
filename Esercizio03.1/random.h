#ifndef __Random__
#define __Random__
#include "function.h"

class Random {

private:
  	int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
  	
public:
	Random();
	~Random();
	void RandGenSetup(); // Added in order to set the initial seed
	void SetRandom(int*, int, int);
	void SaveSeed();
	
	double Uniform(void); // Uniform distribution between 0 and 1
	double Uniform(double min, double max);
	double Linear(double a, double b, double m, double q);
	double Gauss();
	double Gauss(double mean, double sigma); // Box-Muller Method
	double Exp(double lambda);
	double Lorentz(double mean, double gamma);
	double Cos();
	double Angle();
	bool Bool();
	
	// Inverse Sampling and Accept/Reject Methods (Classic, Hybrid)
	double InverseSampling(Function* F);
	double AR(Function* f, double xmin, double xmax, double ymax);
	double AR(Function* f, Function* inv_F);
};

#endif
