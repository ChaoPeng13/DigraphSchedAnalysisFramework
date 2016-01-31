/* \file SchedulabilityAnalysis.h
 * this file implements the schedulability analysis of digraphs and statflows
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 07-sept-2015 : initial revision (CP)
 *
 */

#ifndef SCHEDULABILITYANALYSIS_H_
#define SCHEDULABILITYANALYSIS_H_

//#include <vld.h>

#include "Digraph.h"
#include "Stateflow.h"

class SchedulabilityAnalysis {
public:
	//===================================================================================
	// Digraph information
	//===================================================================================
	// Statistics Data
	static double tfRatio0;
	static double tfRatio1;

	static int nIrreducibleDigraphs;
	// Statistics Time
	static double tDigraphCalSum;
	static double tDigraphCalLinearBounds;
	static double tDigraphCalTF0;
	static double tDigraphCalTF1;
	static double tDigraphCalTF2;
	static double tDigraphRBF;
	static double tDigraphIBF;
	// Vectors
	static vector<double> vec_tfRatio0;
	static vector<double> vec_tfRatio1;

	static vector<int> vec_nIrreducibleDigraphs;

	static vector<double> vec_tDigraphCalSum;
	static vector<double> vec_tDigraphCalLinearBounds;
	static vector<double> vec_tDigraphCalTF0;
	static vector<double> vec_tDigraphCalTF1;
	static vector<double> vec_tDigraphCalTF2;
	static vector<double> vec_tDigraphRBF;
	static vector<double> vec_tDigraphIBF;
	//===================================================================================
	// Stateflow information
	//===================================================================================
	static int nStateflows;
	static double totalUtilization;
	static double avgDegree;
	// statistics Results
	static int nExactStaticOffset;
	static int nExactArbitraryOffset;
	static int nRBFStaticOffset;
	static int nRBFArbitraryOffset;
	static int nRBFArbitraryOffsetBySimpleDigraph;
	static int nRBFArbitraryOffsetByPreciseDigraph;
	static int nIBFStaticOffset;
	static int nIBFArbitraryOffset;
	static int nIBFArbitraryOffsetBySimpleDigraph;
	static int nIBFArbitraryOffsetByPreciseDigraph;
	static int nLinearUpperRBF;
	static int nLinearUpperIBF;
	static int nLinearUpperCSUM;
	// Statistics Data
	static double bpRatio0; // c^rbf+c^dbf /c^sum 
	static double bpRatio1; // c^ibf+c^dbf / c^sum

	static int nIrreducibleStateflows;
	// Statistics Time
	static double tPrepareStateflow; // time for preparing stateflow such as calculating gcd and hyperperiod, generating exeuction request matrix and etc. 

	static double tCalCSum; // time for calculating c^sum
	static double tCalLinearBounds; // time for calculating the linear upper bounds
	static double tCalTF0; // time for calculating t_f baed on c^sum
	static double tCalTF1; // time for calculating t_f based on c^rbf and c^dbf
	static double tCalTF2; // time for calculating t_f based on c^ibf and c^dbf
	static double tCalNonPeriod; // time for calculating the matrix power without periodicity property
	static double tCalPeriod; // time for calculating the matrix power with periodicity property
	static double tCalDiffPeriod; // the different time for calculating the matrix power
	static double tCalRequestExecutionMatrixWithoutPeriodicityProperty;
	static double tCalRequestExecutionMatrixWithPeriodicityProperty;
	static double tGenerateCriticalActionPairs;
	static double tGenerateRequestFunctions;
	static double tGenerateRequestFunctionAbstractTree;
	
	static double tExactStaticOffset;
	static double tExactArbitraryOffset;
	static double tRBFStaticOffset;
	static double tRBFArbitraryOffset;
	static double tRBFArbitraryOffsetBySimpleDigraph;
	static double tRBFArbitraryOffsetByPreciseDigraph;
	static double tIBFStaticOffset;
	static double tIBFArbitraryOffset;
	static double tIBFArbitraryOffsetBySimpleDigraph;
	static double tIBFArbitraryOffsetByPreciseDigraph;
	static double tLinearUpperRBF;
	static double tLinearUpperIBF;
	static double tLinearUpperCSUM;
	
	// Vectors
	static vector<int> vec_nStateflows;
	static vector<double> vec_totalUtilization;
	static vector<double> vec_avgDegree;

	static vector<int> vec_nExactStaticOffset;
	static vector<int> vec_nExactArbitraryOffset;
	static vector<int> vec_nRBFStaticOffset;
	static vector<int> vec_nRBFArbitraryOffset;
	static vector<int> vec_nRBFArbitraryOffsetBySimpleDigraph;
	static vector<int> vec_nRBFArbitraryOffsetByPreciseDigraph;
	static vector<int> vec_nIBFStaticOffset;
	static vector<int> vec_nIBFArbitraryOffset;
	static vector<int> vec_nIBFArbitraryOffsetBySimpleDigraph;
	static vector<int> vec_nIBFArbitraryOffsetByPreciseDigraph;
	static vector<int> vec_nLinearUpperRBF;
	static vector<int> vec_nLinearUpperIBF;
	static vector<int> vec_nLinearUpperCSUM;

	static vector<double> vec_bpRatio0;
	static vector<double> vec_bpRatio1;
	static vector<int> vec_nIrreducibleStateflows;

	static vector<double> vec_tPrepareStateflow;

	static vector<double> vec_tCalCSum;
	static vector<double> vec_tCalLinearBounds;
	static vector<double> vec_tCalTF0;
	static vector<double> vec_tCalTF1;
	static vector<double> vec_tCalTF2;
	static vector<double> vec_tCalNonPeriod;
	static vector<double> vec_tCalPeriod;
	static vector<double> vec_tCalDiffPeriod;
	static vector<double> vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty;
	static vector<double> vec_tCalRequestExecutionMatrixWithPeriodicityProperty;
	static vector<double> vec_tGenerateCriticalActionPairs;
	static vector<double> vec_tGenerateRequestFunctions;
	static vector<double> vec_tGenerateRequestFunctionAbstractTree;
	
	static vector<double> vec_tExactStaticOffset;
	static vector<double> vec_tExactArbitraryOffset;
	static vector<double> vec_tRBFStaticOffset;
	static vector<double> vec_tRBFArbitraryOffset;
	static vector<double> vec_tRBFArbitraryOffsetBySimpleDigraph;
	static vector<double> vec_tRBFArbitraryOffsetByPreciseDigraph;
	static vector<double> vec_tIBFStaticOffset;
	static vector<double> vec_tIBFArbitraryOffset;
	static vector<double> vec_tIBFArbitraryOffsetBySimpleDigraph;
	static vector<double> vec_tIBFArbitraryOffsetByPreciseDigraph;
	static vector<double> vec_tLinearUpperRBF;
	static vector<double> vec_tLinearUpperIBF;
	static vector<double> vec_tLinearUpperCSUM;
                      
	//=====================================================================================
	// Method start
	//=====================================================================================
	static void save_results();
	static void set_zero();
	static void reset();
	static void output_one_vector(ofstream& fout, string sVec, vector<double> vec);
	static void output_one_vector(ofstream& fout, string sVec, vector<int> vec);
	static void output_vectors(ofstream& fout);
	/* Choice used to identify which busy length to choose 
	 * 0 -> original busy length, i.e. csum
	 * 1 -> C^rbf and C^dbf
	 * 2 -> C^ibf and C^dbf
	 */
	
	// prepare t_f 
	static void prepare_all_digraphs(Digraph** digraphs, int n);
	static void generate_critical_vertices(Digraph** digraphs, int n);

	static bool rbf_analysis(Digraph** digraphs, int n, int choice);
	static bool rbf_analysis(Digraph** digraphs, int i, int t, int tprim);

	static bool ibf_analysis(Digraph** digraphs, int n, int choice);
	static bool ibf_analysis(Digraph** digraphs, int i, int t, int tprim);

	/* Choice used to identify which busy length to choose 
	 * 0 -> original busy length, i.e. csum
	 * 1 -> C^rbf and C^dbf
	 * 2 -> C^ibf and C^dbf
	 */

	// prepare linear upper bounds and the maximal time length
	static void prepare_all_stateflows(Stateflow** sfs, int n, int choice);

	// reset some containters for all the stateflows
	static void reset_calculating_containers(Stateflow** sfs, int n);

	// calculate the request execution matrices, in particular to statistic the runtimes
	static void calculate_request_execution_matrices_without_periodicity_property(Stateflow** sfs, int n, int choice);
	static void calculate_request_execution_matrices_with_periodicity_property(Stateflow** sfs, int n, int choice);

	static void generate_critical_action_pair(Stateflow** sfs, int n);

	static bool generate_request_function(Stateflow** sfs, int n, bool output);
	static bool generate_request_function(Stateflow* sf, int i, int maxDeadline, bool output);

	static void intra_loop_dominance(vector<RequestFunction>& vec_rf);
	static void end_loop_dominance(vector<RequestFunction>& vec_rf);

	static void generate_request_function_abstract_tree(Stateflow** sfs, int n, bool output);
	static void generate_request_function_abstract_tree(Stateflow** sf, int i, int stime, int deadline, bool output);
	static void generate_request_function_abstract_tree(Stateflow* sf, int stime, int deadline, bool output);

	static bool isAbstract(vector<RFNode*> vec_rfnode);

	// Exact schedulability analysis
	static bool exact_sched_analysis(Stateflow** sfs, int n, int choice,bool output, ostream& out);
	static bool exact_sched_analysis(Stateflow** sfs, int i, int s, vector<ActionPair> vec_ap);
	static bool exact_sched_analysis(Stateflow** sfs, int i, int s, ActionPair ap);
	static bool exact_sched_analysis(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int stime, int deadline, int gcd);
	static bool exact_sched_analysis(vector<vector<RFNode*>> Omega, int wcet, int stime, int deadline, int gcd);
	static bool exact_sched_analysis(vector<RFNode*> vec_rfnode, int wcet, int stime, int deadline, int gcd);
	static int calculate_remaining_workload(vector<RFNode*> vec_rfnode, int stime, int gcd);

	// Exact analysis with arbitrary offset
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int i, vector<ActionPair> vec_ap);
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int i, ActionPair ap);
	static bool exact_analysis_arbitrary_offset(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int deadline, int gcd);
	static bool exact_analysis_arbitrary_offset(vector<vector<RFNode*>> Omega, int wcet, int deadline, int gcd);
	static bool exact_analysis_arbitrary_offset(vector<RFNode*> vec_rfnode, int wcet, int deadline, int gcd);

	// Approximate schedulability analysis based on RBF with static offsets
	static bool rbf_analysis_static_offset(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_static_offset_index(Stateflow** sfs, int i);
	static bool rbf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);
	static double calculate_remaining_rbf_workload(Stateflow** sfs, int i, int stime, int gcd);

	// Approximate schedulability analysis based on IBF with static offsets
	static bool ibf_analysis_static_offset(Stateflow** sfs, int n, int choice);
	static bool ibf_analysis_static_offset_index(Stateflow** sfs, int i);
	static bool ibf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);
	static double calculate_remaining_ibf_workload(Stateflow** sfs, int i, int stime, int gcd);

	// Approximate schedulability analysis based on RBF with arbitrary offsets
	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i);
	static bool rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline);

	// RBF analysis on simple models
	static bool rbf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_arbitrary_offset_by_simple_digraphs(Digraph** digraphs, int i);
	static bool rbf_analysis_arbitrary_offset_by_simple_digraphs(Digraph** digraphs, int i, int wcet, int deadline);
	// RBF analysis on precise models
	static bool rbf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice);
	
	// Approximate schedulability analysis based on IBF with arbitrary offsets
	static bool ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i);
	static bool ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline);

	// IBF analysis on simple models
	static bool ibf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice);
	// IBF analysis on precise models
	static bool ibf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice);

	// Linear schedulability analysis based on RBF with static offset
	static bool lu_rbf_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_rbf_sched_analysis_index(Stateflow** sfs, int i);
	static bool lu_rbf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);

	// Linear schedulability analysis based on IBF with static offset
	static bool lu_ibf_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_ibf_sched_analysis_index(Stateflow** sfs, int n);
	static bool lu_ibf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);

	// Linear schedulability analysis based on CSUM with static offset
	static bool lu_csum_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_csum_sched_analysis_index(Stateflow** sfs, int n);
	static bool lu_csum_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);
	
	// output
	static void output_critical_action_pair(Stateflow** sfs, int n, ostream& out);

	static void output_critical_request_function(Stateflow** sfs, int n, ostream& out);
	static void output_critical_request_function(Stateflow* sf, ostream& out);
	static void output_critical_request_function(vector<RequestFunction> vec_rf, ostream& out);

	static void output_request_function_abstract_tree(Stateflow** sfs, int n, ostream& out);
	static void output_request_function_abstract_tree(Stateflow* sf, ostream& out);
	static void output_request_function_abstract_tree(AbstractRequestFunctionTree arft, ostream& out);
	
};

#endif