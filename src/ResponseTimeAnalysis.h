/* \file ResponseTimeAnalysis.h
 * this file implements the response time analysis of digraphs and statflows
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 25-sept-2015 : initial revision (CP)
 *
 */

#ifndef RESPONSETIMEANALYSIS_H_
#define RESPONSETIMEANALYSIS_H_

#include "Digraph.h"
#include "Stateflow.h"

class ResponseTimeAnalysis {
public:
	static void calculate_all_response_times(Stateflow** sfs, int n);

	static void calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index);
	static void calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index, Transition* tran);
	static double calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index, Transition* tran, int start, int gcd);

	static void calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index);
	static void calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index, Transition* tran);
	static double calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index, Transition* tran, int s);

	static void calculate_response_time_by_rbf_with_arbitrary_offsets(Stateflow** sfs, int index);
	static void calculate_response_time_by_rbf_with_arbitrary_offsets(Stateflow** sfs, int index, Transition* tran);

	static void calculate_response_time_by_ibf_with_arbitrary_offsets(Stateflow** sfs, int index);
	static void calculate_response_time_by_ibf_with_arbitrary_offsets(Stateflow** sfs, int index, Transition* tran);

};

#endif