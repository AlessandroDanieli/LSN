#ifndef __Position__
#define __Position__
#include <vector>
#include <string>
using namespace std;

class Position {
protected:
	vector<double> X;
	
public:
	Position();
	Position(double, double); // 2D x, y
	Position(double, double, double); // 3D x, y, z
	Position(vector<double> P);
	Position(const Position&);
	~Position();
	Position& operator=(const Position&);
	
	//## Set Methods ##//
	void Set(double , double); // 2D Set method
	void Set(double , double , double ); // 3D Set method
	void Set(vector<double> );
	
	//## Get Methods ##//
	vector<double> GetCoord() const { return X; }
	unsigned int GetDim() const { return X.size(); }
	double GetR() const;
	double GetX() const;
	double GetY() const;
	double GetZ() const;
	double GetRho() const;
	double GetPhi(const string& opt="") const;
	double GetTheta(const string& opt="") const;
	
	//## Other Methods ##//
	vector<double> Angles(const string& opt="") const;// N-dim angular coordinates: x1 = rsin(t1)...sin(tN), x2 = rsin(t1)...cos(tN), x3 = rsin(t1)...cos(t(N-1)) ... xN = rcos(t1)
	double Distance() const;
	double Distance(const Position&) const;
	void Print() const;
	void PrintPolar(const string& opt="") const; // 2D polar coordinates
	void PrintCyl(const string& opt="") const;
	void PrintSph(const string& opt="") const;
};
#endif
