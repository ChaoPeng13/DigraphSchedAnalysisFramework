/* \file GraphAlgorithms.cpp
*  this file implements a series of graph algorithms. 
*  \author Chao Peng
*  
*  Changes
*  ------
*  25-aug-2015 : initial revision (CP)
*
*/

#ifndef GRAPHALGORITHMS_H_
#define GRAPHALGORITHMS_H_

//#include <vld.h>

#include "Digraph.h"
#include "GeneralDirectedGraph.h"
#include "Stateflow.h"

using namespace std;

/// In this class, all the member functions are static for supporting direct access.
class GraphAlgorithms {
public:

	static double calculate_maximum_cycle_mean(Digraph* digraph);
	/*
	 * \brief An implementation of the iterative algorithm for calculating utilization
	 * in Stigge et al. The digraph real-time task model. RTAS2011.
	 * In case of the large number of nodes in the unit digraph.
	 */
	static double calculate_maximum_cycle_mean2(Digraph* digraph);
	static double calculate_maximum_cycle_mean(double** A, int n);
	static double calculate_maximum_cycle_mean(GeneralDirectedGraph* gDigraph);
	

	/// /brief generate strongly connected componets for the digraph
    /// The algorithm is implemented after "Cormen et al: Introduction to
	/// agorithms", Chapter 22.5. It has a running time of O(V + E)
	static vector<Digraph*> generate_strongly_connected_components(Digraph* digraph);

	/// /brief an implementation of the depth-first-search algorithm after "Cormen et al: Introduction to
	/// agorithms", Chapter 22.3. The vector new_node_vec is used to specify the searching order, which 
	/// should be used in the step 3 of the strongly connected component algorithm.
	static void depth_first_search(Digraph* digraph, vector<Node*> new_node_vec);
	static void depth_first_search_visit(Node* node, int &time);

	/// Depth-first-search on the tranpose graph
	static void reverse_depth_first_search(Digraph* digraph, Node** node_array, int n, vector<Digraph*>& sccs);
	static void reverse_depth_first_search_visit(Node* node, Digraph* digraph);

	// /brief generate strongly connected componets for the digraph
    /// The algorithm is implemented after "Cormen et al: Introduction to
	/// agorithms", Chapter 22.5. It has a running time of O(V + E)
	static vector<GeneralDirectedGraph*> generate_strongly_connected_components(GeneralDirectedGraph* gDigraph);

	/// /brief an implementation of the depth-first-search algorithm after "Cormen et al: Introduction to
	/// agorithms", Chapter 22.3. The vector new_node_vec is used to specify the searching order, which 
	/// should be used in the step 3 of the strongly connected component algorithm.
	static void depth_first_search(GeneralDirectedGraph* gDigraph, vector<Node*> new_node_vec);
	//static void depth_first_search_visit(Node* node, int &time);

	/// Depth-first-search on the tranpose graph
	static void reverse_depth_first_search(GeneralDirectedGraph* gDigraph, Node** node_array, int n, vector<GeneralDirectedGraph*>& sccs);
	static void reverse_depth_first_search_visit(Node* node, GeneralDirectedGraph* digraph);

	/// Generate execution request matrix for a digraph
	static double** generate_exec_request_matrix(Digraph* digraph);
	static double** generate_exec_request_matrix(Digraph* digraph,set<int> iSet);

	/// Calculate the sum of wcets of all the nodes in the digraph
	static void calculate_csum(Digraph* digraph);

	/// Calculate the linear bounds derived from rbf and ibf
	static void calculate_tight_linear_bounds(Digraph* digraph);

	/// Generate a simple digraph from the stateflow
	static Digraph* generate_simple_digraph(Stateflow* sf);

	/// Generate a precise digraph from the stateflow
	static Digraph* generate_precise_digraph(Stateflow* sf);

	/// An implementation of http://www.geeksforgeeks.org/detect-cycle-in-a-graph/
	/// Detect cycle in a digraph, DFS
	static bool isCyclic(Digraph* digraph);
	static bool isCyclicUtil(Digraph* digraph, int v, bool* visited, bool* recStack);

	/// Detect cycle in a stateflow, DFS
	static bool isCyclic(Stateflow* sf);
	static bool isCyclicUtil(Stateflow* sf, int v, bool* visited, bool* recStack);
};

#endif