/* \file StateAndTransition.h
*  this file decribes the State and Transition class. 
*  Stateflows (Finite State Machines) consist of a number of states and transitions
*  \author Chao Peng
*  
*  Changes
*  ------
*  31-Aug-2015 : initial revision (CP)
*
*/

#ifndef STATEANDTRANSITION_H_
#define STATEANDTRANSITION_H_

//#include <vld.h>

#include <string>
#include <vector>
#include "Utility.h"

using namespace std;

class State;
class Transition;

class State {
public:
	string name;
	int index;
	list<Transition*> in;
	list<Transition*> out;

	bool visited; // used to detect cycle

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	State(string _name, int _index) {
		name = _name;
		index = _index;
	}
	
	~State() {
		in.clear();
		out.clear();
	}

	std::string toString() {
		return name+"-"+Utility::int_to_string(index);
	}
};

class Transition {
public:
	State* src;
	State* snk;
	int period;
	int priority;
	int wcet;
	int scale_wcet; // might be used to do breakdown factors or action extensibility
	int deadline;
	string name;

	bool visited;

	// Response times
	map<int,double> responseTimeByRbfWithStaticOffsets;
	map<int,double> responseTimeByIbfWithStaticOffsets;

	double responseTimeByRbfWithArbitraryOffsets;
	double responseTimeByIbfWithArbitraryOffsets;

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	Transition(State* _src, State* _snk) : src(_src), snk(_snk) {}
	Transition(State* _src, State* _snk, string _name) : src(_src), snk(_snk), name(_name) {}

	string toString() {
		if (name.empty())
			return src->toString()+"->"+snk->toString();
		else 
			return src->toString()+"->"+snk->toString()+"\t"+name;
	}

	// Sorted by the priority
	bool operator < (const Transition& rhs) const {
		return priority < rhs.priority;
	}
};

#endif