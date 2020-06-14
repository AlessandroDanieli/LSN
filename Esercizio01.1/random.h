#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void RandGenSetup(); // Added in order to set the initial seed
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Uniform(void); // Uniform distribution between 0 and 1
  double Uniform(double min, double max);
  double Linear(double a, double b, double m, double q);
  double Gauss(double mean, double sigma); // Box-Muller Method
  double Exp(double lambda);
  double Lorentz(double mean, double gamma);
};

#endif // __Random__
