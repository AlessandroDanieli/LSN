#ifndef __Vector__
#define __Vector__
#include <vector>
#include <string>
using namespace std;

class Vector {
private:
	double* V = nullptr;
	double norm = -1;
	unsigned int d = 0;
	
public:
	Vector();
	Vector(unsigned int _d);
	Vector(const Vector& x);
	Vector(const vector<double> x);
	//Vector(double radius, const vector<double> angles = vector<double>()); // Angular coordinates
	~Vector();
	
	//## Set Methods ##//
	void FillValue(unsigned int n, double value); // Clears the vector and fills a new one with n copies of value
	void FillValue(double value); // Sets each entry = value (size unchanged)
	void SetComponent(unsigned int i, double x);
	void Set(Vector& x);
	void Set(double radius, vector<double> angles);
	void Append(double value);
	
	//## Get Methods ##//
	unsigned int GetDim() const { return d; }
	double* GetV() const { return V; }
	double GetComponent(unsigned int i) const;
	double GetThetaSph();
	double GetPhiSph();
	
	//## Copy Methods ##//
	void Copy(vector<double> x);
	void Copy(Vector);
		
	//## Manipulation Methods ##//
	vector<double> Angles(const string& opt="");
	double Norm();
	double Distance(const Vector&) const;
	double Distance2(const Vector&) const;
	void Normalize();
	
	double operator[](unsigned int i) const;
	double& operator[](unsigned int i);
	bool operator==(const Vector& X) const;
	Vector& operator=(const Vector& X);
	Vector operator+(const Vector& Y) const;
	void operator-();
	Vector operator-(const Vector& Y) const;
	void operator+=(const Vector& Y);
	void operator-=(const Vector& Y);
	double operator*(const Vector& Y) const;
	
	Vector operator*(const double scalar);
	void operator*=(const double scalar);
	Vector operator/(const double scalar);
	void operator/=(const double scalar);
	
	Vector Parallel(const Vector& Y) const;
	Vector Normal(const Vector& Y) const;
	
	void Swap(unsigned int p, unsigned int q);
	void Print(unsigned int prec = 5) const;
	void PrintSph(const string& opt="") ;
	void Clear();
};

#endif

