#ifndef REQUESTFUNCTION_H_
#define REQUESTFUNCTION_H_

//#include "StateAndTransition.h"
#include "Stateflow.h"
#include "LiftingPoint.h"

class RequestFunction {
public:
	int priority;
	int gcd;
	int hyperperiod;
	int start; // start time of rf
	int finish; // finish time of rf, [s,f)

	int lastValue; // the last value
	Transition* lastTran; // the last transition
	bool isCounted; // whether the last transition has been counted 

	list<LiftingPoint> lps; // lifting points
	list<Transition*> trans; // list of transitions
	map<int,map<int,int>> record; // rf[s,f)
	map<int,int> rf; // rf(t)

	RequestFunction() {}
	~RequestFunction() {
		lps.clear();
	}

	RequestFunction(int _priority, int _gcd, int _hyperperiod, int _start, int _finish) {
		priority = _priority;
		gcd = _gcd;
		hyperperiod = _hyperperiod;
		start = _start;
		finish = _finish;

		lastValue = 0;
		lastTran = NULL;
		
		LiftingPoint lp = LiftingPoint(start,lastValue);
		lps.push_back(lp);
	}

	RequestFunction(const RequestFunction& rf) {
		priority = rf.priority;
		gcd = rf.gcd;
		hyperperiod = rf.hyperperiod;
		start = rf.start;
		finish = rf.finish;

		lastValue = rf.lastValue;
		lastTran = rf.lastTran;
		isCounted = rf.isCounted;

		lps = rf.lps;
		trans = rf.trans;
	}

	void update_next_action (int curTime, Transition* t) {
		trans.push_back(t);
		if (lastTran == NULL) {
			lastTran = t;
			isCounted = false;

			if (curTime != start) {
				cerr << "rf start time error." <<endl;
				exit(EXIT_FAILURE);
			}
			return;
		}
		
		finish = curTime;

		if (isCounted) {
			LiftingPoint& lp = lps.back();
			lp.x = curTime;
		}
		else {
			lastValue += lastTran->wcet;
			LiftingPoint lp = LiftingPoint(curTime, lastValue);
			lps.push_back(lp);
			isCounted = true;
		}

		if (t != NULL) {
			lastTran = t;
			isCounted = false;
		}
	}

	int findRecord(int s, int f) {
		if (record.find(s) == record.end())
			return -1;
		if (record[s].find(f) == record[s].end())
			return -1;
		return record[s][f];
	}

	int request_function(int s, int f) {
		int interval = f-s;
		s = s%hyperperiod;
		f = s+interval;
		if ( f <= start) return 0;

		int ret = findRecord(s,f);
		if (ret != -1) return ret;

		int sValue = -1;
		for (list<LiftingPoint>::iterator iter = lps.begin(); iter != lps.end(); iter++) {
			LiftingPoint lp = *iter;
			if (lp.x >= s) {
				sValue = lp.y;
				break;
			}
		}
		if (sValue == -1) return 0;

		int fValue = -1;
		for (list<LiftingPoint>::iterator iter = lps.begin(); iter != lps.end(); iter++) {
			LiftingPoint lp = *iter;
			if (lp.x >= f) {
				fValue = lp.y;
				break;
			}
		}

		if (fValue == -1) fValue = lastValue;

		record[s][f] = fValue-sValue;
		return fValue-sValue;
	}

	int request_function(int t) {
		t = ceil(1.0*t/gcd)*gcd;
		
		if (rf.find(t) != rf.end()) return rf[t];

		int value = INT_MIN;
		for (int s=0; s<=finish; s+=gcd) {
			int f = s+t;
			value = max(value, request_function(s,f));
		}

		rf[t] = value;
		return value;
	}

	int get_length() {
		return finish - start;
	}

	const bool operator>(RequestFunction& rf) {
		if (priority != rf.priority) return false;
		if (start != rf.start) return false;

		for (int t = start; t<=finish; t+=gcd)
			if (request_function(start,t) < rf.request_function(start,t)) return false;

		return true;
	}

	void output(ostream& out) {
		out << "Request Function: <" << priority << "," << gcd << ","
			<< start << "," << finish << "," << lastValue << ">" <<endl;

		out << "Lifting Ponts:" <<endl;
		for (list<LiftingPoint>::iterator iter = lps.begin(); iter != lps.end(); iter++) {
			LiftingPoint lp = *iter;
			out << "=>";
			lp.output(out);
		}
		out << endl;

		out << "Transition List" << endl;
		out << trans.size() <<endl;
		for (list<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter ++) {
			Transition* t = *iter;
			out << "=>";
			if (t == NULL) out <<t;
			else {
				out << t->toString();
			}
			out << endl;
		}
		out << endl;
	}

};

#endif