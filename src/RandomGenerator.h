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
	static enum PERIOD_CHOICE
	{
		ExactAnalysis0, // {10, 20, 25, 50} and {1,2}
		ExactAnalysis1, // {5, 10, 20, 25, 50} and {1,2,4}
		ExactAnalysis2, // {5, 10, 20, 25, 50} and {1,2}
		ExactAnalysis3, // {10, 20, 25, 50} and {1,2,4}
		ApproximateAnalysis0, // {5, 10, 20, 25, 50, 100}
		ApproximateAnalysis1, // {5, 10, 20, 25, 50, 100} and {1,2,4,5}
		ApproximateAnalysis2, // {5, 10, 20, 25, 50, 100} and {1,2,4}
	};

	static int STATEFLOW_PERIOD[6]; // {5, 10, 20, 25, 50, 100};
	static int STATEFLOW_PERIOD2[5]; // {10, 20, 25, 50, 100}
	static int STATEFLOW_PERIOD3[4]; // {10, 20, 25, 50}
	static int STATEFLOW_PERIOD4[5]; // {5, 10, 20, 25, 50}
	//static int STATEFLOW_FACTOR[6]; // {1, 2, 4, 5, 10, 20}
	static int STATEFLOW_BASE[8]; // {1100, 1300,2300, 3100, 3700, 3000, 5000, 7000}
	//static int STATEFLOW_BASE[5]; // {1100, 1300,2300, 3100, 3700, 3000, 5000}
	static int STATEFLOW_FACTOR[5]; // {1, 2, 4, 5, 10, 20}
	static int STATEFLOW_FACTOR2[4]; // {1,2,4,5}
	static int STATEFLOW_FACTOR3[3]; // {1,2,4}
	static int STATEFLOW_FACTOR4[2]; // {1,2}

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
	/// Generate one random stateflow with nState states and nTran transitions
	/// and make sure the out-degree of any state <= 2
	static Stateflow* generate_one_stateflow2(int index, int scale, int nState, int nTran);

	static void setup_periods_for_one_stateflow_for_exact_analysis_0(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_1(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_2(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_3(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_0(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_1(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_2(Stateflow* sf, int scale);

	static void setup_wcets_for_one_stateflow(Stateflow* sf, double util);

	static Stateflow* generate_one_stateflow_with_util(int index, int scale, double util, int nState, int nTran);
	static Stateflow* generate_one_stateflow_with_util2(int index, int scale, double util, int nState, int nTran);
	/// Generate one random stateflow system with num stateflows and tUtil (total) utilization
	static Stateflow** generate_stateflow_system(int num, int maxState, double tUtil);

	/// Generate one random stateflow system based on STATEFLOW_FACTOR and STATEFLOW_BASE
	static Stateflow** generate_stateflow_system2(int num, int maxState, double tUtil);

	/// Generate one random stateflow system (two partitions: 1~num-1 and num)
	static Stateflow** generate_stateflow_system3(int num, int maxState, double tUtil);

	/// Generate one random stateflow system (the last one is special)
	static Stateflow** generate_stateflow_system4(int num, int maxState, double tUtil);
	static Stateflow** generate_stateflow_system5(int num, int maxState, double tUtil);
	/// if period_choice == 1, we use the function setup_periods_for_one_stateflow
	/// else if period_choice == 2, we use the function setup_periods_for_one_stateflow2
	/// else if period_choice == 3, we use the function setup_periods_for_one_stateflow3
	static Stateflow** generate_stateflow_system6(int num, int maxState, double tUtil,int period_choice);
	/// if period_choice == 1, we use the function setup_periods_for_one_stateflow
	/// else if period_choice == 2, we use the function setup_periods_for_one_stateflow2
	/// else if period_choice == 3, we use the function setup_periods_for_one_stateflow3
	static Stateflow** generate_stateflow_system7(int num, int maxState, double tUtil, int period_choice);

	/// Corresponding to the generate_stateflow_system6, we add one parameter scc_probability
	static Stateflow** generate_stateflow_system8(int num, int maxState, double tUtil,int period_choice, double scc_probability);

	static Stateflow** generate_stateflow_system_for_exact_analysis(int num, int maxState, double tUtil,int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability);

	/// Calculate the number of edges such taht the average degree for a vertex is around 3
	static int calculate_num_edge(int numNode);
};

#endif