#include "position.h"
#include "needle.h"
#include <cmath>
using namespace std;

Needle :: Needle(): Position() {}
Needle :: ~Needle() {}

Needle :: Needle(double x, double y): Position(x, y) {}

Needle :: Needle(double x, double y, double _l, double _angle): Position(x, y) {
	l = _l;
	angle = _angle;
}

void Needle :: Set(double _x, double _y, double _l, double _angle) {
	X[0] = _x;
	X[1] = _y;
	l = _l;
	angle = _angle;
	return;
}

bool Needle :: IsCrossed(double min, double max, double d) {
	double xmax = X[0]+l/2.*abs(cos(angle));
	double xmin = X[0]-l/2.*abs(cos(angle));
	unsigned int n = (unsigned int)(max-min)/d;
	for (unsigned int i = 0; i <= n; i++) {
		if (min+i*d >= xmin && min+i*d <= xmax) return true;
	}
	return false;
}

bool Needle :: IsCrossed(vector<double>& x_lines) {
	double dx = l/2.*abs(cos(angle));
	double xmax = X[0]+dx;
	double xmin = X[0]-dx;
	for (double x : x_lines) {
		if (x >= xmin && x <= xmax) return true;
	}
	return false;
}

