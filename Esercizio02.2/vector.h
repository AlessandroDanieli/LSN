#ifndef __Vector__
#define __Vector__
#include <vector>
#include <string>
using namespace std;

class Vector {
private:
	vector<double> V;
	
public:
	Vector();
	Vector(const Vector& x);
	Vector(const vector<double> x);
	Vector(double radius, const vector<double> angles = vector<double>()); // Angular coordinates
	~Vector();
	
	//## Set Methods ##//
	void FillValue(unsigned int n, double value); // Clears the vector and fills a new one with n copies of value
	void FillValue(double value); // Sets each entry = value (size unchanged)
	void SetComponent(unsigned int i, double x);
	void Set(Vector x);
	void Set(vector<double> x);
	void Set(double radius, vector<double> angles);
	void Append(double value);
	
	//## Get Methods ##//
	unsigned int GetDim() const { return V.size(); }
	vector<double> GetV() const { return V; }
	double GetComponent(unsigned int i) const;
	
	//## Copy Methods ##//
	void Copy(vector<double>&);
	void Copy(Vector);
		
	//## Manipulation Methods ##//
	vector<double> Angles(const string& opt="") const;
	double Norm(const Vector&) const;
	double Norm() const;
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
	
	Vector operator*(const double scalar) const;
	void operator*=(const double scalar);
	Vector operator/(const double scalar) const;
	void operator/=(const double scalar);
	
	Vector Parallel(const Vector& Y) const;
	Vector Normal(const Vector& Y) const;
	
	void Swap(unsigned int p, unsigned int q);
	void Print() const;
	void PrintSph(const string& opt="") const;
	void Clear() { V.clear(); }
};

#endif

