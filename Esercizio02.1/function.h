#ifndef __Function__
#define __Function__
#include <cmath>

class Function {
public:
	virtual double Eval(double x) const = 0;
};

// Generic function (integrals)
class Integrand : public Function {
private:
public:
	Integrand() {};
	~Integrand() {};
	double Eval(double x) const override { return M_PI/2.*cos(M_PI/2.*x); }
	double operator()(double x) const { return Eval(x); }
};

class Linear : public Function {
protected:
	double m, q;
public:
	Linear(): m(0), q(1) {};
	Linear(double _m, double _q): m(_m), q(_q) {};
	~Linear() {};
	void Set(double _m, double _q) { m = _m; q = _q; }
	double GetM() const { return m; }
	double GetQ() const { return q; }
	double Eval(double x) const { return m*x + q; }
	double operator()(double x) const { return Eval(x); }
};

// Linear Importance Sampling correctly normalized 
class LinearIS : public Linear {
private:
	double min, max; // Used in Importance Sampling methods
public:
	LinearIS(): Linear(), min(0), max(1) {};
	LinearIS(double _m, double _q, double _min, double _max): Linear(_m, _q), min(_min), max(_max) {};
	~LinearIS() {};
	void Set(double _m, double _q, double _min, double _max) { m = _m; q = _q; min = _min; max = _max; }
	double GetM() const { return m; }
	double GetQ() const { return q; }
	double GetMin() const { return min; }
	double GetMax() const { return max; }
	double Eval(double x) const override { return (m*x+q)/((max-min)*(m*(max+min)/2.+q)); }
	double operator()(double x) const { return Eval(x); }
};

class Gauss : public Function {
private:
	double mean, sigma, k;
public:
	Gauss(): mean(0), sigma(1), k(1) {}; // Multiplicative coefficient
	Gauss(double m, double s); // Normalized distribution
	Gauss(double m, double s, double _k): mean(m), sigma(s), k(_k) {};
	~Gauss() {};
	
	void SetMean(double m) { mean = m; }
	void SetSigma(double s) { sigma = s; }
	void SetK(double _k) { k = _k; }
	double GetMean() const { return mean; }
	double GetSigma() const { return sigma; }
	double GetK() const { return k; }
	void Normalize();
	double Eval(double x) const override; 
	double operator()(double x) const { return Eval(x); }
};

class Lorentz : public Function {
private:
	double mean, gamma, k;
public:
	Lorentz() {};
	Lorentz(double m, double s); // Normalized distribution
	Lorentz(double m, double g, double _k): mean(m), gamma(g), k(_k) {};
	~Lorentz() {};
	
	void SetMean(double m) { mean = m; }
	void SetSigma(double g) { gamma = g; }
	void SetK(double _k) { k = _k; }
	double GetMean() const { return mean; }
	double GetGamma() const { return gamma; }
	double GetK() const { return k; }
	void Normalize();
	double Eval(double x) const override { return k / ( (x-mean)*(x-mean) + gamma*gamma ); }
	double operator()(double x) const { return Eval(x); }
};

#endif
