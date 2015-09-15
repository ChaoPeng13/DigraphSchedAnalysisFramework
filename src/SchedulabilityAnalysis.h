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
	// Statistics Data
	static double bpRatio0; // c^rbf+c^dbf /c^sum 
	static double bpRatio1; // c^ibf+c^dbf / c^sum
	// Statistics Time
	static double tCalCSum; // time for calculating c^sum
	static double tCalLinearBounds; // time for calculating the linear upper bounds
	static double tCalTF0; // time for calculating t_f baed on c^sum
	static double tCalTF1; // time for calculating t_f based on c^rbf and c^dbf
	static double tCalTF2; // time for calculating t_f based on c^ibf and c^dbf
	static double tRBFStaticOffset;  // time for RBF satic offset
	static double tRBFArbitraryOffset; // time for RBF arbitrary offset
	static double tIBFStaticOffset; // time for IBF satic offset
	static double tIBFArbitraryOffset; // time for IBF arbitrary offset

	//=====================================================================================
	// Method start
	//=====================================================================================

	static bool rbf_analysis(Digraph** digraphs, int n);
	static int get_rbf_value(Digraph** digraphs, int i, int t, int tprim);

	static bool ibf_analysis(Digraph** digraphs, int n);
	static int get_ibf_value(Digraph** digraphs, int i, int t, int tprim);

	/* Choice used to identify which busy length to choose 
	 * 0 -> original busy length, i.e. csum
	 * 1 -> C^rbf and C^dbf
	 * 2 -> C^ibf and C^dbf
	 */

	static bool rbf_analysis_static_offset(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_static_offset(Stateflow** sfs, int i, int s, int t, int tprim);

	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim);

	static bool ibf_analysis_static_offset(Stateflow** sfs, int n, int choice);
	static double get_ibf_static_offset(Stateflow** sfs, int i, int s, int t, int tprim);

	static bool ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static double  get_ibf_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim);
};

#endif