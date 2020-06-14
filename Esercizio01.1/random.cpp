#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

// Saves the seeds used: allows to continue a calculation starting from the last seed used (not necessarily the original seed) //
void Random :: SaveSeed() {
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Uniform(double min, double max){
   return min+(max-min)*Uniform();
}

double Random :: Linear(double a, double b, double m, double q) {
	if (a == b || m*a+q < 0 || m*b +q < 0) {
		cerr << "[Random :: Linear] Error: Invalid input" << endl;
		return 0;
	}
	else if (m == 0) {
		return Uniform(a, b);
	}
	else {
		double u = Uniform(a, b);
		double r = q/m;
		double t = b+a+2.*r;
		return -r - sqrt(r*r + t*u-a*b);
	}
}

double Random :: Gauss(double mean, double sigma) { // ## Box-Muller Method ## //
   double s=Uniform();
   double t=Uniform();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double lambda) {
	double s=Uniform();
	return -log(s)/lambda; // Inverse of the cumulative distribution
}

double Random :: Lorentz(double mean, double gamma) {
	double s=Uniform();
	return gamma*tan(M_PI*(s-0.5))+mean;
}


// Uniform distribution between 0 and 1
double Random :: Uniform(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

// Added || General setup for the random generator (extracted from main)
void Random :: RandGenSetup() {
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()) {
	  Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	  while ( !input.eof() ) {
		 input >> property;
		 if( property == "RANDOMSEED" ){
		    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		    SetRandom(seed,p1,p2); // ## Modify this if you want to export it to another file: rnd.SetRandom(..) ## //
		 }
	  }
	  input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	return;
}

