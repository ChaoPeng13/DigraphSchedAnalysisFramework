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

#include "Digraph.h"
#include "Stateflow.h"

class SchedulabilityAnalysis {
public:
	static bool rbf_analysis(Digraph** digraphs, int n);
	static int get_rbf_value(Digraph** digraphs, int i, int t, int tprim);

	static bool ibf_analysis(Digraph** digraphs, int n);
	static int get_ibf_value(Digraph** digraphs, int i, int t, int tprim);

	static bool rbf_analysis_static_offset(Stateflow** sfs, int n);
	static bool rbf_analysis_static_offset(Stateflow** sfs, int i, int s, int t, int tprim);

	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int n);
	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim);

	static bool ibf_analysis_static_offset(Stateflow** sfs, int n);
	static int get_ibf_static_offset(Stateflow** sfs, int i, int s, int t, int tprim);

	static bool ibf_analysis_arbitrary_offset(Stateflow** sfs, int n);
	static int  get_ibf_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim);
};

#endif