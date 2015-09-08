/* \file Stateflow.h
 *  This file decribes the Stateflow class.
 *  \author Chao Peng
 *  
 *  Changes
 *  ------
 *  xx-Aug-2015 : Initial revision (CP)
 *
 */

#ifndef STATEFLOW_H_
#define STATEFLOW_H_

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>

#include "Digraph.h"
#include "GeneralDirectedGraph.h"
#include "Utility.h"
#include "StateAndTransition.h"

#pragma once;

using namespace std;

class StartFinishTime {
public:
	int start;
	int finish;
	double value;

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	StartFinishTime(int s, int f, double v) {
		start = s;
		finish = f;
		value = v;
	}

	std::string toString() {
		return Utility::int_to_string(start)+","+Utility::int_to_string(finish);
	}
};

/// \brief Stateflow or Finite State Machine (FSM) class
/// The detailed definition of Stateflows can be found in
/// Zeng and Di Natale, Schedulability analysis of periodic tasks implementing synchronous finite state machines, ECRTS2012
/// Here, we consider of simple (flat) FSMs, without concurrency and hierarchy, as well as other Stateflow extensions and notational conveniences.
/// It should be extended to more general FSMs.
class Stateflow {
public:
	//friend class Digraph;
	int index;

	vector<State*> states;
	vector<Transition*> trans;
	map<State*, int> state_index;
	map<int, State*> index_state;
	map<Transition*, int> tran_index;
	map<int, Transition*> index_tran;

	int iState; // global variable for indexing new state
	int iTran;  // global varibale for indexing new transition

	int scale;

	int gcd; // the greatest common divisor of the periods of the transitions
	int t_gcd; // the greastest common divisor of the periods and wcets of the transisitions 
	int hyperperiod; // the least common multiple of the periods of the transitions
	
	int n_state;

	set<int> rbf_time_instances; // time instances for rbf
	map<int,int> index_time; // time index, first: index, second: time
	map<int,int> time_index; // reverse time index, first: time, second: index
	int n_time_instance; // size of rbf_time_instances
	double**** rbfs; // four-dimension matrix. rbf[i][j][s][f] = rbf_{i,j}[s,f)


	set<int> ibf_time_instances; // time instances for ibf
	double**** ibfs; // four-dimension matrix. ibf[i][j][s][f] = ibf_{i,j}[s,f)

	map<StartFinishTime*, double> rbf_hyperperiod; // rbf of s and f within one hyperperiod  
	map<StartFinishTime*, double> ibf_hyperperiod; // ibf of s and f within one hyperperiod

	double** exec_req_matrix; // execution request matrix x^{k}_{i,j}=rbf_{i,j}[0,kH_F)
	map<int, double**> exec_req_matrix_power; // execution request matrix

	int gper; // generalized period
	int gdef; // generalized defect
	double** gfac; // generalized factor, if reducible, then all the elements are linear factor
	bool isIrred; // reducible of the execution request matrix
	double lfac; // linear factor

	double csum; // linear bound, used to calculate busy-period length
	double tf0; // by csum
	double max_tf0;

	Digraph* simple_digraph;
	Digraph* precise_digraph;
	GeneralDirectedGraph* exec_digraph; // digraph for execution request matrix

	// tighter linear bound
	double crbf;
	double cibf;
	double cdbf;

	double tf1; // by crbf and cdbf
	double max_tf1; 

	double tf2; // by cibf and cdbf
	double max_tf2;

	vector<StartFinishTime*> rbf_vec; // vector of rbf[s,f)

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	Stateflow(int _scale) { 
		this->scale = _scale;
		iState = 0;
		iTran = 0;

		isIrred = false;
	}

	Stateflow(int _index, int _scale) {
		index = _index;
		scale = _scale;
		iState = 0;
		iTran = 0;

		isIrred = false;
	}

	~Stateflow(){
		/*
		state_index.clear();
		index_state.clear();
		tran_index.clear();
		index_tran.clear();
		// release states
		for (vector<State*>::iterator iter = states.begin(); iter != states.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		states.clear();

		// release transitions
		for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		trans.clear();

		rbf_time_instances.clear();
		index_time.clear();
		time_index.clear();

		// delete rbfs
		for (int i=0; i<n_state; i++) {
			for (int j=0; j<n_state; j++) {
				for (int s=0; s< n_time_instance; s++) {
					delete rbfs[i][j][s];
				}
				delete rbfs[i][j];
			}
			delete rbfs[i];
		}
		delete rbfs;

		rbf_hyperperiod.clear();
		
		// delete execution request matrix
		for(map<int, double**>::iterator iter = exec_req_matrix_power.begin(); iter != exec_req_matrix_power.end(); iter++) {
			for (int i=0; i<n_state; i++) {
				delete (iter->second)[i];
			}
			delete iter->second;
		}
		exec_req_matrix_power.clear();

		delete simple_digraph;
		delete precise_digraph;
		delete exec_digraph;
		*/

	}

	void add_state(State* state); // add new state
	void add_transition(Transition* tran); // add new transition

	bool less_compare_trans(const Transition* t1, const Transition* t2);

	void prepare(); // prepare everything

	void calculate_gcd(); // calculate the greatest common divisor of periods of transitions
	void calculate_t_gcd(); // calculate the greatest common divisor of periods and wcets of transitions
	void calculate_hyperperiod(); // calculate the hyperperiod of the stateflow, i.e., the least common multiple of periods of transitions

	void set_state_number(); // set n_state

	void generate_rbf_time_instances();
	void generate_ibf_time_instances();

	void generate_rbfs();
	void generate_ibfs();

	void generate_exec_req_matrix();
	void calculate_exec_req_matrix_power(int tf);

	// Calculate generialized periodicity parameters
	void calculate_generialized_factor();
	void calculate_generialized_period();
	void calculate_generialized_defect();

	// Check whether the execution request matrix is irreducible
	void check_irreducible();

	// Calculate the linear factor
	void calculate_linear_factor();

	// rescale wcets
	void scale_wcet(double factor);

	// Calculate csum and busyperiod length
	void calculate_csum();
	void calculate_tf0(Stateflow** stateflows, int i);

	// Generate a simple digraph for stateflow
	void generate_simple_digraph();
	void generate_precise_digraph();

	/// Calculate the tigher linear upper bounds
	/// true: operating on precise digraph
	/// false: operating on simple digraph
	void calculate_linear_upper_bounds(bool precise);
	void calculate_tf1(Stateflow** stateflows, int i);
	void calculate_tf2(Stateflow** stateflows, int i);

	/// Static offset
	double get_rbf(int start, int finish); // return rbf[s,f)
	double calculate_rbf_within_one_hyperperiod(int start, int finish);
	double calculate_rbf_within_one_hyperperiod(int i, int j, int start, int finish);
	double calculate_rbf_within_multiple_hyperperiods(int start, int finish);
	

	double get_ibf(int start, int finish); // return ibf[s,f)
	double calculate_ibf_within_one_hyperperiod(int start, int finish);
	double calculate_ibf_within_one_hyperperiod(int i, int j, int start, int finish);
	double calculate_ibf_within_multiple_hyperperiods(int start, int finish);

	double get_dbf(int start, int finish); // return dbf[s,f)

	// Arbitrary offset
	double get_rbf(int t); // return rbf(t)
	double get_ibf(int t); // return ibf(t)
	double get_dbf(int t); // return dbf(t)

	/// \brief write a stateflow object into an output stream in graphviz dot format
	void write_graphviz(ostream& out);
};

#endif