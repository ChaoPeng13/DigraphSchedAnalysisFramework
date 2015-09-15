/* \file RandomGenerator.cpp
*  this file implements a random generator for Digraph or Stateflow systemss. 
*  \author Chao Peng
*  
*  Changes
*  ------
*  04-Sept-2015 : initial revision (CP)
*
*/

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

//#include <vld.h>

#include "Digraph.h"
#include "Stateflow.h"

using namespace std;

class RandomGenerator {
public:
	/// Graph properties
	static int DIGRAPH_SCALE;
	static int DIGRAPH_PERIOD[3][2]; // {{2,4},{6,12},{5,10}};
	static int DIGRAPH_SCALE_FACTOR[5]; // {1,2,4,5,10};

	/// Stateflow properties
	static int STATEFLOW_SCALE;
	static int STATEFLOW_PERIOD[6]; // {5, 10, 20, 25, 50, 100};

	/// Generate one random digraph with nNode nodes and nEdge edges
	/// which is the same as the paper ECRTS2013
	static Digraph* generate_one_digraph(int index, int scale, int nNode, int nEdge);
	static Digraph* generate_one_digraph(int index, int scale, int nNode);
	/// Generate one random digraph systems with num digraphs and tUtil (total) utilization
	static Digraph** generate_digraph_system(int num, int maxNode, double tUtil);

	static int calculate_base();
	static int calculate_separation_time(int base);

	/// Generate one random stateflow with nState states and nTran transitions
	/// which is the same as the paper ECRTS2012
	static Stateflow* generate_one_stateflow(int index, int scale, int nState, int nTran);
	/// Generate one random stateflow systems with num stateflows and tUtil (total) utilization
	static Stateflow** generate_stateflow_system(int num, int maxState, double tUtil);

	/// Calculate the number of edges such taht the average degree for a vertex is around 3
	static int calculate_num_edge(int numNode);
};

#endif