#include "ResponseTimeAnalysis.h"

void ResponseTimeAnalysis::calculate_all_response_times(Stateflow** sfs, int n) {
	for (int i=0; i<n; i++) {
		calculate_response_time_by_rbf_with_static_offsets(sfs,i);
		//calculate_response_time_by_ibf_with_static_offsets(sfs,i);
		calculate_response_time_by_rbf_with_arbitrary_offsets(sfs,i);
		//calculate_response_time_by_ibf_with_arbitrary_offsets(sfs,i);
	}
}

void ResponseTimeAnalysis::calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index) {
	Stateflow* sfi = sfs[index];

	for (vector<Transition*>::iterator iter = sfi->trans.begin(); iter != sfi->trans.end(); iter++)
		calculate_response_time_by_rbf_with_static_offsets(sfs,index,*iter);
}

void ResponseTimeAnalysis::calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index, Transition* tran) {
	if (index == 0) {
		tran->responseTimeByRbfWithStaticOffsets.insert(pair<int,double>(0,(double)tran->wcet));
		return;
	}

	// calculate hyperperiod for 0~index stateflows
	int hyperperiod = 1;
	int gcd = 0;
	for (int j=0; j<=index; j++) {
		Stateflow* sfj = sfs[j];
		hyperperiod = Utility::math_lcm(hyperperiod, sfj->hyperperiod);
		gcd = Utility::math_gcd(gcd, sfj->gcd);
	}

	// calculate the set of start times within one hyperperiod
	set<int> sSet;
	for (int s=0; s<=hyperperiod; s+=tran->period)
		sSet.insert(s);

	for (set<int>::iterator iter = sSet.begin(); iter != sSet.end(); iter++) {
		int s = *iter;
		tran->responseTimeByRbfWithStaticOffsets[s] = calculate_response_time_by_rbf_with_static_offsets(sfs, index, tran, s, gcd);
	}
}

double ResponseTimeAnalysis::calculate_response_time_by_rbf_with_static_offsets(Stateflow** sfs, int index, Transition* tran, int s, int gcd) {
	// calculate response times
	int total = 0;
	int f = s;
	while (true) {
		total = tran->wcet;
		for (int j=0; j<index; j++) {
			Stateflow* sfj = sfs[j];
			double rbfj = sfj->get_rbf(s,f);
			total += rbfj;
		}
		if (total <= f-s) break;
		f += gcd;
	}
	
	return total;
}

void ResponseTimeAnalysis::calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index) {
	Stateflow* sfi = sfs[index];

	for (vector<Transition*>::iterator iter = sfi->trans.begin(); iter != sfi->trans.end(); iter++)
		calculate_response_time_by_ibf_with_static_offsets(sfs,index,*iter);
}

void ResponseTimeAnalysis::calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index, Transition* tran) {
	if (index == 0) {
		tran->responseTimeByIbfWithStaticOffsets[0] = tran->wcet;
		return;
	}

	// calculate hyperperiod for 0~index stateflows
	int hyperperiod = 1;
	for (int j=0; j<=index; j++) {
		Stateflow* sfj = sfs[j];
		hyperperiod = Utility::math_lcm(hyperperiod, sfj->hyperperiod);
	}

	// calculate the set of start times within one hyperperiod
	set<int> sSet;
	for (int s=0; s<=hyperperiod; s+=tran->period)
		sSet.insert(s);

	for (set<int>::iterator iter = sSet.begin(); iter != sSet.end(); iter++) {
		int s = *iter;
		tran->responseTimeByIbfWithStaticOffsets[s] = calculate_response_time_by_ibf_with_static_offsets(sfs, index, tran, s);
	}
}

double ResponseTimeAnalysis::calculate_response_time_by_ibf_with_static_offsets(Stateflow** sfs, int index, Transition* tran, int s) {
	// calculate response times
	int f = s+tran->wcet;
	double total = tran->wcet;
	
	while (true) {
		for (int j=0; j<index; j++) {
			Stateflow* sfj = sfs[j];
			double ibfj = sfj->get_ibf(s,f);
			total += ibfj;
		}
		if (total <= f-s)
			break;
		if (f == s+ceil(total)) f++;
		else f = s+ceil(total);
	}
	return total;
}

void ResponseTimeAnalysis::calculate_response_time_by_rbf_with_arbitrary_offsets(Stateflow** sfs, int index) {
	Stateflow* sfi = sfs[index];

	for (vector<Transition*>::iterator iter = sfi->trans.begin(); iter != sfi->trans.end(); iter++)
		calculate_response_time_by_rbf_with_arbitrary_offsets(sfs,index,*iter);
}

void ResponseTimeAnalysis::calculate_response_time_by_rbf_with_arbitrary_offsets(Stateflow** sfs, int index, Transition* tran) {
	if (index == 0) {
		tran->responseTimeByRbfWithArbitraryOffsets = tran->wcet;
		return;
	}

	// calculate gcd for 0~index stateflows
	int gcd = 0;
	for (int j=0; j<=index; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd, sfj->gcd);
	}

	// calculate response times
	int total = 0;
	int t = gcd;
	while (true) {
		total = tran->wcet;
		for (int j=0; j<index; j++) {
			Stateflow* sfj = sfs[j];
			double rbfj = sfj->get_rbf(t);
			total += rbfj;
		}
		if (total <= t) break;
		t += gcd;
	}
	
	tran->responseTimeByRbfWithArbitraryOffsets = total;
}

void ResponseTimeAnalysis::calculate_response_time_by_ibf_with_arbitrary_offsets(Stateflow** sfs, int index) {
	Stateflow* sfi = sfs[index];

	for (vector<Transition*>::iterator iter = sfi->trans.begin(); iter != sfi->trans.end(); iter++)
		calculate_response_time_by_ibf_with_arbitrary_offsets(sfs,index,*iter);
}

void ResponseTimeAnalysis::calculate_response_time_by_ibf_with_arbitrary_offsets(Stateflow** sfs, int index, Transition* tran) {
	if (index == 0) {
		tran->responseTimeByIbfWithArbitraryOffsets = tran->wcet;
		return;
	}

	int t = tran->wcet;
	double total = tran->wcet;
	
	while (true) {
		for (int j=0; j<index; j++) {
			Stateflow* sfj = sfs[j];
			double ibfj = sfj->get_ibf(t);
			total += ibfj;
		}
		if (total <= t)
			break;
		if (t == ceil(total)) t++;
		else t = ceil(total);
	}
	
	tran->responseTimeByIbfWithArbitraryOffsets = total;
}