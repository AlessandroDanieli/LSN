#ifndef __Needle__
#define __Needle__
#include "position.h"
using namespace std;

class Needle : public Position {
private:
	double l, angle; // length, length x projection
public:
	Needle();
	~Needle();
	Needle(double x, double y);
	Needle(double x, double y, double _l, double _angle);
	
	void Set(double _x, double _y, double _l, double _angle);

	double GetLength() { return l; }
	double GetAngle() { return angle; }
	Position GetPosition() const { return X; }
	
	bool IsCrossed(vector<double>& line_x_pos);
	bool IsCrossed(double min, double max, double d);
};

#endif
