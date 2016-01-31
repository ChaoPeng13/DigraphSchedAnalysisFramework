#ifndef LITFTINGPOINT_H_
#define LITFTINGPOINT_H_

class LiftingPoint {
public:
	int x;
	int y;

	LiftingPoint(int _x, int _y) {
		x = _x;
		y = _y;
	}
	~LiftingPoint() {
		//cout << "destruct LiftingPoint" <<endl;
	}

	bool operator>(const LiftingPoint &lp) {
		if (x <= lp.x && y >= lp.y) return true;
		return false;
	}

	bool operator==(const LiftingPoint &lp) {
		if ( x==lp.x && y==lp.y) return true;
		return false;
	}

	void output(ostream& out) {
		out << "<" << x << "," << y << ">"<<endl;
	}
};

#endif

