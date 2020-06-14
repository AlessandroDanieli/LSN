#ifndef __Randomwalk__
#define __Randomwalk__
#include "random.h"
#include "vector.h"
	
using namespace std;

class Randomwalk {
protected:
	Random* Rnd;
	Vector P;

public:
	Randomwalk();
	Randomwalk(Random* r);
	Randomwalk(Random* r, vector<double> position);
	Randomwalk(Random* r, Vector position);
	~Randomwalk();
	
	void SetRandomGen(Random* r) { Rnd = r; }
	void SetPosition(Vector position) { P = position; }
	
	Random* GetRandomGen() const { return Rnd; }
	Vector GetPosition() const { return P; }
	
	void Move(Vector v) { P += v; }
	void RandomMove(double distance, const unsigned int k = 1); // Move k times
	double Distance() const;
	double Distance(const vector<double>&);
	double Distance(const Vector&);
	double Distance(const Randomwalk&);
	
	void PrintPosition() const;
};


class RandomwalkLattice : public Randomwalk {
private:
	vector<Vector> Gen;
	
public:
	RandomwalkLattice();
	RandomwalkLattice(Random* r, Vector position, vector<Vector> generators);
	~RandomwalkLattice();
	
	void SetGenerators(const vector<Vector>& generators);
	void RandomMove(const unsigned int k = 1);
};

#endif
