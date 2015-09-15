/* \file Digraph.h
*  this file decribes the digraph class. 
*  directed graphs consist of a number of nodes and edges
*  \author Chao Peng
*  
*  Changes
*  ------
*  21-Aug-2015 : initial revision (CP)
*
*/

#ifndef DIGRAPH_H_
#define DIGRAPH_H_

//#include <vld.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <list>

#include "Utility.h"
#include "NodeAndEdge.h"

#pragma once

extern double POS_INFINITY;
extern double NEG_INFINITY;

using namespace std;

class UnitDigraph;
class GranDigraph;

/// \brief an implementation of the real-time task model. 
/// M. stigge et al., the digraph real-time task model, rtas2011
class Digraph {
public:
	int index;
	int scale; // scale the (small) wcet of any node to be an integer

	// now it is enough to decribe a directed graph by nodes and edges.
	vector<Node*> node_vec;
	vector<Edge*> edge_vec;

	// identify nodes and edges
	map<Node*, int> node_to_index;
	map<int, Node*> index_to_node;
	map<Edge*, int> edge_to_index;
	map<int, Edge*> index_to_edge;

	int iNode;
	int iEdge;

	// gcd
	int gcd; // greatest common divisor of the minimum separation times 
	vector<Digraph*> sccs; // the vector of all the strongly connected components
	vector<Digraph*> hccs; // the vector of all the highly connected components

	/// we have to explicity define two transformed digraph: UDRT and GDRT.
	UnitDigraph* unit_digraph;
	GranDigraph* gran_digraph;
	
	// linear periodicity parameters
	// i guess these parameters respectively derived from rbf and ibf might be same.
	double linear_factor;

	// linear upper bound parameters
	double c_rbf;
	double c_ibf;
	double c_dbf;
	double c_sum;

	// linear upper bounds
	double tf0; // by c_sum
	double tf1;  // zeng and di natale, using max-plus algebra to improve the analysis of non-cyclic task models, ecrts2013
	double tf2;  // new linear bound

	int tf; // the maximal instance for calculated rbf and ibf to do schedulability analysis

	// properties of the directed graph
	bool strongly_connected;
	int maximum_degree;
	double average_degree;

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	Digraph() {
		index = 0;
		scale = 0;
		unit_digraph = NULL;
		gran_digraph = NULL;

		iNode = 0;
		iEdge = 0;
	}

	Digraph(int _index) {
		index = _index;
		scale = 0;
		unit_digraph = NULL;
		gran_digraph = NULL;

		iNode = 0;
		iEdge = 0;
	}

	Digraph(int _index, int _scale) { //default constructor
		index = _index;
		scale = _scale;
		unit_digraph = NULL;
		gran_digraph = NULL;

		iNode = 0;
		iEdge = 0;

		// todo: we assert that node_vec and edge_vec should be empty here.
	}
	~Digraph();

	void add_node(Node* p_node) { 
		this->node_vec.push_back(p_node); 
		node_to_index[p_node] = iNode; 
		index_to_node[iNode] = p_node;
		iNode++;
	}

	void add_edge(Edge* p_edge) { 
		this->edge_vec.push_back(p_edge); 
		edge_to_index[p_edge] = iEdge;
		index_to_edge[iEdge] = p_edge;
		p_edge->src_node->out.push_back(p_edge);
		p_edge->snk_node->in.push_back(p_edge);
		iEdge++;
	}

	void add_scc_edge(Edge* p_edge) {
		this->edge_vec.push_back(p_edge); 
		edge_to_index[p_edge] = iEdge;
		index_to_edge[iEdge] = p_edge;
		p_edge->src_node->scc_out.push_back(p_edge);
		p_edge->snk_node->scc_in.push_back(p_edge);
		iEdge++;
	}

	void calculate_gcd() {
		int _gcd = 0;
		for (vector<Edge*>::iterator iter = this->edge_vec.begin(); iter != this->edge_vec.end(); iter++) {
			Edge* edge = *iter;
			_gcd = Utility::math_gcd(_gcd, edge->separationTime); 
		}
		this->gcd = _gcd;
	}

	void prepare_digraph();

	void generate_strongly_connected_components();
	void generate_highly_connected_components();
	void check_strongly_connected();

	void calculate_linear_factor();
	void calculate_linear_upper_bounds();

	void scale_wcet(double factor);

	void calculate_tf0(Digraph** digraphs, int i);
	void calculate_tf1(Digraph** digraphs, int i);
	void calculate_tf2(Digraph** digraphs, int i);
	

	void prepare_rbf_calculation(bool debug);
	void prepare_ibf_calculation(bool debug);

	/// \brief return the rbf(t)
	double rbf(int t);
	/// \brief return the rbf execution request matrix at t
	double** rbf_exec_req_matrix(int t, int& n);

	/// \brief return the ibf(t)
	double ibf(int t);
	/// \brief return the ibf execution request matrix at t
	double** ibf_exec_req_matrix(int t, int& n);

	/// \brief return the dbf(t)
	double dbf(int t);

	/// \brief write a digraph object into an output stream in graphviz dot format
	void write_graphviz(ostream& out);
};

/// \brief an implmentation of the unit digraph task. 
/// Zeng and Di Natale, Using Max-Plus Algebra to Improve the Analysis of Non-cyclic Task Models, ECRTS2013
class UnitDigraph: public Digraph {
public:
	// point to the original digraph
	Digraph* origin;

	int n_size; // the size of nodes
	int scale;
	int gcd;  // greatest common divisor of the minimum separation times

	set<int> iSet; // set of indices for node->wcet != 0
	map<Node*, int> original_node_to_index; // the index of the original node
	map<int, Node*> index_to_original_node; // the original node at the index

	/// the execution request matrix of the transformed unit digraph task (udrt)
	/// with the greatest common divisor of the minimum separation times 
	/// (or the periods of all the edges) but not 1.
	double** matrix; // two-dimension matrix
	std::map<int, double**> matrix_map; // matrix power
	std::map<int, int> maximum_element_map; // the maximum element of the execution request matrix at the special time

	double lfac;
	int lper;
	int ldef;

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	UnitDigraph(Digraph* digraph) {
		origin = digraph;
		scale = digraph->scale;
		lfac = digraph->linear_factor;
		gcd = digraph->gcd;
		tf = digraph->tf;

		unit_digraph = NULL;
		gran_digraph = NULL;

		iNode = 0;
		iEdge = 0;

		// TODO: we assert that node_vec and edge_vec should be empty here.
	}
	~UnitDigraph();

	void add_node(Node* p_node) { 
		this->node_vec.push_back(p_node); 
		node_to_index[p_node] = iNode; 
		index_to_node[iNode] = p_node;
		if (p_node->wcet != 0) iSet.insert(iNode);
		iNode++;
	}

	void add_node(Node* p_node, Node* o_node) { 
		this->node_vec.push_back(p_node); 
		node_to_index[p_node] = iNode; 
		index_to_node[iNode] = p_node;
		if (p_node->wcet != 0) {
			iSet.insert(iNode);
			original_node_to_index[o_node] = iNode;
			index_to_original_node[iNode] = o_node;
		}
		iNode++;
	}

	void add_edge(Edge* p_edge) { 
		this->edge_vec.push_back(p_edge); 
		edge_to_index[p_edge] = iEdge;
		index_to_edge[iEdge] = p_edge;
		iEdge++;
	}

	void calculate_gcd() {
		int _gcd = 0;
		for (vector<Edge*>::iterator iter = origin->edge_vec.begin(); iter != origin->edge_vec.end(); iter++) {
			Edge* edge = *iter;
			_gcd = Utility::math_gcd(_gcd, edge->separationTime); 
		}
		this->gcd = _gcd;
	}

	void prepare(bool debug);
	void generate_unit_digraph();
	void scc_generate_unit_digraph(); // generate unit digraph from a strongly connected component
	void generate_exec_request_matrix();
	void calculate_linear_period(bool debug);
	void calculate_linear_defect();

	/// \brief return the rbf execution request matrix at the time t power
	void calculate_exec_request_matrix_power(int tf);
	int get_rbf(int t);

	/// \brief calculate the dbf_{i,j}(t) = dbf(v_i,v_j,t) by the l-MAD property assumption
	/// dbf(v_i,v_j,t) = rbf(v_i,vj,t-d(v_j)-\min\limits_{v_k:(v_k,v_j)\in\mathbb{E}} p(v_k,v_j)) + e(v_j)
	int calculate_demand_bound_function(int i, int j, int t);
	int get_dbf(int t);
};

/// \brief an implmentation of the granularity digraph task. 
/// RTNS2015 version
class GranDigraph: public Digraph {
public:
	// point to the original digraph
	Digraph* origin;

	int n_size; // the size of nodes
	int scale;
	int gcd;  // greatest common divisor of the minimum separation times and wcets

	/// the execution request matrix of the transformed granularity digraph task (udrt)
	/// with the greatest common divisor of the minimum separation times 
	/// (or the periods of all the edges) and wcets but not 1.
	double** matrix; // two-dimension matrix
	std::map<int, double**> matrix_map; // matrix power
	std::map<int, int> maximum_element_map; // the maximum element of the execution request matrix at the special time

	double lfac;
	int lper;
	int ldef;

	///=============================================================================================================================
	/// method functions
	///=============================================================================================================================
	GranDigraph(Digraph* digraph) {
		origin = digraph;
		scale = digraph->scale;
		lfac = digraph->linear_factor;
		tf = digraph->tf;

		unit_digraph = NULL;
		gran_digraph = NULL;

		iNode = 0;
		iEdge = 0;

		// TODO: we assert that node_vec and edge_vec should be empty here.
	}
	~GranDigraph();

	void add_node(Node* p_node) { this->node_vec.push_back(p_node); node_to_index.insert(pair<Node*,int>(p_node,iNode++));}
	void add_edge(Edge* p_edge) { this->edge_vec.push_back(p_edge); edge_to_index.insert(pair<Edge*,int>(p_edge,iEdge++));}

	void calculate_gcd() {
		// calculate the greatest common divisor of separation times and wcets
		int _gcd = 0;
		for (vector<Edge*>::iterator iter = origin->edge_vec.begin(); iter != origin->edge_vec.end(); iter++) {
			Edge* edge = *iter;
			_gcd = Utility::math_gcd(_gcd, edge->separationTime); 
		}

		for (vector<Node*>::iterator iter = origin->node_vec.begin(); iter != origin->node_vec.end(); iter++) {
			Node* node = *iter;
			_gcd = Utility::math_gcd(_gcd, node->wcet);
		}

		this->gcd = _gcd;
	}

	void prepare(bool debug);
	void generate_gran_digraph();
	void generate_exec_request_matrix();
	void calculate_linear_period(bool debug);
	void calculate_linear_defect();

	/// /brief return the ibf execution request matrix at the time t power
	void calculate_exec_request_matrix_power(int tf);
	int get_ibf(int t);
};

#endif