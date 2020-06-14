#ifndef __Brownian__
#define __Brownian__
#include "random.h"
using namespace std;

class Brownian { // Time independent parameters <mean, sigma>
protected:
	Random* rnd;
	double w, t; // Starting parameters
	double mean, sigma;

public:
	Brownian();
	Brownian(Random* r);
	Brownian(Random* r, double w0, double t0, double _mean, double _sigma);
	~Brownian();
	
	void SetRandGen(Random* r);
	void SetW(double w0) { w = w0; }
	void SetT(double t0) { t = t0; }
	void SetMean(double m) { mean = m; }
	void SetSigma(double s) { sigma = s; }
	
	Random* GetRandGen() const { return rnd; }
	double GetW() const { return w; }
	double GetT() const { return t; }
	double GetMean() const { return mean; }
	double GetSigma() const { return sigma; }
	
	void Move(); // Assuming mean = 0, variance = Delta(t) = 1
	void Move(double delta_t); // variance = Delta(t) = t_<i> - t_<i-1>
};


class GBM : public Brownian {
private:
	double S;
public:
	GBM();
	GBM(Random* r);
	GBM(Random* r, double mean, double sigma, double S0);
	~GBM();
	
	void SetS(double s) { S = s; }
	void Set(double _t, double _S) { t = _t; S = _S; }
	double GetS() const { return S; }
	
	void Evolve(double delta_t, unsigned int steps = 1);
};

#endif
