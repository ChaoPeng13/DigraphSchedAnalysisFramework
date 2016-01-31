#ifndef ACTIONPAIR_H_
#define ACTIONPAIR_H_

#include "StateAndTransition.h"

class ActionPair {
public:
	int e; // wcet
	int d; // deadline
	int stime; // trigger time
	int priority; 
	Transition* t;

	ActionPair(Transition* _t, int _stime) {
		t = _t;
		e = t->wcet;
		stime = _stime;
	}

	ActionPair(int _e, int _d, int _stime, int _priority, Transition* _t) {
		e = _e;
		d = _d;
		stime = _stime;
		priority = _priority;
		t = _t;
	}

	const bool operator>(const ActionPair& ap) {
		if (stime == ap.stime && e >= ap.e && d <= ap.d && priority >= ap.priority) return true;
		return false;
	}

	void output(ostream& out) {
		out << "Transition: " << t->toString() << "\t" << "<" << priority << "," << stime << "," << e << "," << d << ">" <<endl;
	}

	string toString() {
		return "Transition: " + t->toString() + "\t" + "<" + Utility::int_to_string(priority) + "," + Utility::int_to_string(stime) 
			+ "," + Utility::int_to_string(e) + "," + Utility::int_to_string(d) + ">";
	}
};

#endif