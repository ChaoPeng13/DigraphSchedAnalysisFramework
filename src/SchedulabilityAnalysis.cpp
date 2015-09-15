#include "SchedulabilityAnalysis.h"
#include "Timer.h"

#include <algorithm>

double SchedulabilityAnalysis::bpRatio0 = 0;
double SchedulabilityAnalysis::bpRatio1 = 0;

double SchedulabilityAnalysis::tCalCSum = 0;
double SchedulabilityAnalysis::tCalLinearBounds = 0;
double SchedulabilityAnalysis::tCalTF0 = 0;
double SchedulabilityAnalysis::tCalTF1 = 0;
double SchedulabilityAnalysis::tCalTF2 = 0;

double SchedulabilityAnalysis::tRBFStaticOffset = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffset = 0;
double SchedulabilityAnalysis::tIBFStaticOffset = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffset = 0;

bool SchedulabilityAnalysis::rbf_analysis(Digraph** digraphs, int n) {
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		// digraph->calculate_linear_factor();
		digraph->calculate_linear_upper_bounds();
		if (i!=0) {
			digraph->calculate_tf0(digraphs, i);
			digraph->calculate_tf1(digraphs, i);
			digraph->calculate_tf2(digraphs, i);
		} 
		else {
			digraph->tf0 = 0;
			digraph->tf1 = 0;
			digraph->tf2 = 0;
		}
		cout<<"tf0="<<digraph->tf0<<"\ttf1="<<digraph->tf1<<"\ttf2="<<digraph->tf2<<endl;
	}

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int tf2 = ceil(digraph->tf2);
		tf_max = max(tf_max, tf2);
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max/digraph->gcd;
		digraph->prepare_rbf_calculation(false);
	}

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		int tf2 = digraph->tf2;
		if (tf2 == 0) { schedulable[i] = true; continue; }

		for (int t=digraph->gcd; t<tf2; t+=digraph->gcd) {
			schedulable[i] = false;

			int tprim = digraph->gcd;
			while(tprim <= t) {
				int rbf = get_rbf_value(digraphs, i, t, tprim);
				if (rbf <= tprim) {
					schedulable[i] = true;
					break;
				}
				if (tprim == rbf) 
					tprim = rbf+digraph->gcd;
				else tprim = rbf;
			}

			if (!schedulable[i])
				return false;
		}
	}
	
	return true;
}

int SchedulabilityAnalysis::get_rbf_value(Digraph** digraphs, int i, int t, int tprim) {
	double ret = digraphs[i]->dbf(t);

	for (int k=0; k<i; k++) {
		ret += digraphs[k]->rbf(tprim);
	}

	return ceil(ret);
}

bool SchedulabilityAnalysis::ibf_analysis(Digraph** digraphs, int n) {
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		// digraph->calculate_linear_factor();
		digraph->calculate_linear_upper_bounds();
		if (i!=0) {
			digraph->calculate_tf0(digraphs, i);
			digraph->calculate_tf1(digraphs, i);
			digraph->calculate_tf2(digraphs, i);
		} 
		else {
			digraph->tf0 = 0;
			digraph->tf1 = 0;
			digraph->tf2 = 0;
		}
		cout<<"tf0="<<digraph->tf0<<"\ttf1="<<digraph->tf1<<"\ttf2="<<digraph->tf2<<endl;
	}

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int tf2 = ceil(digraph->tf2);
		tf_max = max(tf_max, tf2);
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max/digraph->gcd;
		digraph->prepare_rbf_calculation(false);
		digraph->tf = tf_max/digraph->gran_digraph->gcd;
		digraph->prepare_ibf_calculation(false);
	}

	

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=1; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int gcd = digraph->gran_digraph->gcd;

		int tf2 = digraph->tf2;

		for (int t=gcd; t<tf2; t+=gcd) {
			schedulable[i] = false;

			int tprim = gcd;
			while(tprim <= t) {
				int ibf = get_ibf_value(digraphs, i, t, tprim);
				if (ibf <= tprim) {
					schedulable[i] = true;
					break;
				}
				if (tprim == ibf) 
					tprim = ibf+gcd;
				else tprim = ibf;
			}

			if (!schedulable[i])
				return false;
		}
	}
	
	return true;
}

int SchedulabilityAnalysis::get_ibf_value(Digraph** digraphs, int i, int t, int tprim) {
	double ret = digraphs[i]->dbf(t);

	for (int k=0; k<i; k++) {
		ret += digraphs[k]->ibf(tprim);
	}

	return ceil(ret);
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_t_gcd();

		timer.start();
		sf->calculate_csum();
		timer.end();
		tCalCSum += timer.getTime();

		//sf->write_graphviz(cout);
		//sf->precise_digraph->write_graphviz(cout);

		//sf->check_irreducible();

		sf->calculate_generialized_period();
		sf->calculate_generialized_defect();

		sf->generate_simple_digraph();

		timer.start();
		sf->generate_precise_digraph();
		sf->calculate_linear_upper_bounds(true);
		timer.end();
		tCalLinearBounds += timer.getTime();

		timer.start();
		sf->calculate_tf0(sfs,i);
		timer.end();
		tCalTF0 += timer.getTime();

		timer.start();
		sf->calculate_tf1(sfs,i);
		timer.end();
		tCalTF1 += timer.getTime();

		timer.start();
		sf->calculate_tf2(sfs,i);
		timer.end();
		tCalTF2 += timer.getTime();

		bpRatio0 += sf->tf1/sf->tf0;
		bpRatio1 += sf->tf2/sf->tf0;
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}

	int hyperperiod = 1;
	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		hyperperiod = Utility::math_lcm(hyperperiod, sf->hyperperiod);

		set<int> sSet;
		sSet.insert(0);

		// get the set of start times within one hyperperiod
		for (int j=0; j<=i; j++) {
			Stateflow* hsf = sfs[j];
			/*
			for (int k=1; k<hsf->n_time_instance; k++) {
				int t = hsf->index_time[k];
				for (int l=1; l<hyperperiod/t; l++) sSet.insert(l*t);
			}
			*/
			for (int k=0; k<hyperperiod/hsf->hyperperiod; k++) {
				for (int l=1; l<hsf->n_time_instance; l++) {
					int s = k*hsf->hyperperiod+hsf->index_time[l];
					if (s != hyperperiod)
						sSet.insert(s);
				}
			}
		}

		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf==0) {schedulable[i] = true; continue;}

		for (set<int>::iterator iter = sSet.begin(); iter != sSet.end(); iter++) {
			int s = *iter;
			int f = s + tf;

			// get the set of times for t: s<=t<=f, t is an integer multiple of a deadline
			set<int> tSet;
			for (int deadline = sf->gcd; deadline <=f; deadline+=sf->gcd) {
				if (deadline <=s) continue;
				tSet.insert(deadline);
			}
			tSet.insert(ceil(1.0*f/sf->gcd)*sf->gcd);
			/*
			for (int k=0; k<=f/sf->hyperperiod; k++) {
				for (int j=0; j<sf->n_time_instance; j++) {
					int deadline = k*sf->hyperperiod + sf->index_time[j] + sf->gcd;
					if ( deadline <= s) continue;
					if ( deadline > f) continue;

					tSet.insert(deadline);
				}
			}
			*/
			for (set<int>::iterator iter2 = tSet.begin(); iter2 != tSet.end(); iter2++) {
				int t = *iter2;

				schedulable[i] = false;
				// get the set of times for tprim<=t, tprim is an integer multiple of an instant in s
				set<int> tprimSet;
				for (int j=0; j<=i; j++) {
					Stateflow* hsf = sfs[j];
					for (int tprim = hsf->gcd; tprim<=t; tprim+=hsf->gcd) {
						if (tprim >= s)
							tprimSet.insert(tprim);
					}
				}

				for (set<int>::iterator iter3 = tprimSet.begin(); iter3 != tprimSet.end(); iter3++) {
					int tprim = *iter3;
					if (rbf_analysis_static_offset(sfs, i, s, t, tprim)) {
						schedulable[i] = true;
						break;
					}
				}

				if (!schedulable[i]) {
					delete[] schedulable;
					timer.end();
					tRBFStaticOffset += timer.getTime();
					return false;
				}
			}
		}
	}

	delete[] schedulable;
	timer.end();
	tRBFStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset(Stateflow** sfs, int i, int s, int t, int tprim) {
	double total = 0;
	double dbfi = sfs[i]->get_dbf(s,t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double rbfj = sfs[j]->get_rbf(s,tprim);
		total += rbfj;

		if (total > tprim-s) return false;
	}
	if (total > tprim-s) return false;
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_t_gcd();
		sf->calculate_csum();
		sf->generate_simple_digraph();
		sf->generate_precise_digraph();
		sf->calculate_linear_upper_bounds(true);

		//sf->check_irreducible();

		sf->calculate_generialized_period();
		sf->calculate_generialized_defect();

		sf->calculate_tf0(sfs,i);
		sf->calculate_tf1(sfs,i);
		sf->calculate_tf2(sfs,i);
			
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf == 0) {schedulable[i] = true; continue;}

		set<int> tSet;
		for(int j=0; j<=i; j++) {
			Stateflow* hsf = sfs[j];
			for (int t=hsf->gcd; t<=tf; t+=hsf->gcd)
				tSet.insert(t);
			tSet.insert(hsf->gcd);
		}
		tSet.insert(sf->gcd);

		// get the set of times for t: s<=t<=f, t is an integer multiple of a deadline
		for (set<int>::iterator iter = tSet.begin(); iter != tSet.end(); iter++) {
			int t = *iter;
			schedulable[i] = false;

			set<int> tprimSet;
			for (int j=0; j<=i; j++) {
				Stateflow* hsf = sfs[j];
				for (int tprim = hsf->gcd; tprim<=t; tprim+=hsf->gcd) {
					tprimSet.insert(tprim);
				}
			}

			for (set<int>::iterator iter = tprimSet.begin(); iter != tprimSet.end(); iter++) {
				int tprim = *iter;
				if(rbf_analysis_arbitrary_offset(sfs,i,t,tprim)) {
					schedulable[i] = true;
					break;
				}
			}

			if (!schedulable[i]) {
				delete[] schedulable;
				timer.end();
				tRBFArbitraryOffset += timer.getTime();
				return false;
			}
		}
	}
	delete[] schedulable;
	timer.end();
	tRBFArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim) {
	double total = 0;
	double dbfi = sfs[i]->get_dbf(t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double rbfj = sfs[j]->get_rbf(tprim);
		total += rbfj;

		if (total > tprim) return false;
	}
	if (total > tprim) return false;
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset(Stateflow** sfs, int n, int choice) {
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_t_gcd();
		sf->calculate_csum();
		sf->generate_simple_digraph();
		sf->generate_precise_digraph();
		sf->calculate_linear_upper_bounds(true);

		//sf->check_irreducible();

		sf->calculate_generialized_period();
		sf->calculate_generialized_defect();

		sf->calculate_tf0(sfs,i);
		sf->calculate_tf1(sfs,i);
		sf->calculate_tf2(sfs,i);
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	Timer timer;
	timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}

	int hyperperiod = 1;
	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		hyperperiod = Utility::math_lcm(hyperperiod, sf->hyperperiod);

		set<int> sSet;
		sSet.insert(0);

		// get the set of start times within one hyperperiod
		for (int j=0; j<=i; j++) {
			Stateflow* hsf = sfs[j];
			/*
			for (int k=1; k<hsf->n_time_instance; k++) {
				int t = hsf->index_time[k];
				for (int l=1; l<hyperperiod/t; l++) sSet.insert(l*t);
			}
			*/
			for (int k=0; k<hyperperiod/hsf->hyperperiod; k++) {
				for (int l=1; l<hsf->n_time_instance; l++) {
					int s = k*hsf->hyperperiod+hsf->index_time[l];
					if (s != hyperperiod)
						sSet.insert(s);
				}
			}
		}

		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf==0) {schedulable[i] = true; continue;}

		for (set<int>::iterator iter = sSet.begin(); iter != sSet.end(); iter++) {
			int s = *iter;
			int f = s + tf;

			// get the set of times for t: s<=t<=f, t is an integer multiple of a deadline
			set<int> tSet;
			for (int deadline = sf->gcd; deadline <=f; deadline+=sf->gcd) {
				if (deadline <=s) continue;
				tSet.insert(deadline);
			}
			tSet.insert(ceil(1.0*f/sf->gcd)*sf->gcd);

			for (set<int>::iterator iter2 = tSet.begin(); iter2 != tSet.end(); iter2++) {
				int t = *iter2;

				schedulable[i] = false;
				
				int tprim = s;
				while (tprim <= t) {
					int total = ceil(get_ibf_static_offset(sfs, i, s, t, tprim))+s;
					if (total <= tprim) {
						schedulable[i] = true;
						break;
					}
					if (tprim == total) total++;
					tprim = total;
				}

				if (!schedulable[i]) {
					delete[] schedulable;
					timer.end();
					tIBFStaticOffset += timer.getTime();
					return false;
				}
	
			}
		}
	}
	delete[] schedulable;
	timer.end();
	tIBFStaticOffset += timer.getTime();
	return true;
}

double SchedulabilityAnalysis::get_ibf_static_offset(Stateflow** sfs, int i, int s, int t, int tprim) {
	double total = 0;
	double dbfi = sfs[i]->get_dbf(s,t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double ibfj = sfs[j]->get_ibf(s,tprim);
		total += ibfj;

		//if (total > tprim-s) return false;
	}

	return total;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_t_gcd();
		sf->calculate_csum();
		sf->generate_simple_digraph();
		sf->generate_precise_digraph();
		sf->calculate_linear_upper_bounds(true);

		//sf->check_irreducible();

		sf->calculate_generialized_period();
		sf->calculate_generialized_defect();

		sf->calculate_tf0(sfs,i);
		sf->calculate_tf1(sfs,i);
		sf->calculate_tf2(sfs,i);
			
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	Timer timer;
	timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		int tf;
		if (choice == 0) tf = ceil(sf->tf0);
		else if (choice == 1) tf = ceil(sf->tf1);
		else if (choice == 2) tf = ceil(sf->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf == 0) {schedulable[i] = true; continue;}

		set<int> tSet;
		for(int j=0; j<=i; j++) {
			Stateflow* hsf = sfs[j];
			for (int t=hsf->gcd; t<=tf; t+=hsf->gcd)
				tSet.insert(t);
			tSet.insert(hsf->gcd);
		}
		tSet.insert(sf->gcd);

		// get the set of times for t: s<=t<=f, t is an integer multiple of a deadline
		for (set<int>::iterator iter = tSet.begin(); iter != tSet.end(); iter++) {
			int t = *iter;
			schedulable[i] = false;

			int tprim = 0;
			while (tprim <= t) {
				int ibf = get_ibf_arbitrary_offset(sfs, i, t, tprim);
				if (ibf <= tprim) {
					schedulable[i] = true;
					break;
				}
				if (tprim == ibf) 
					tprim = ibf++;
				else tprim = ibf;
			}

			if (!schedulable[i]) {
				delete[] schedulable;
				timer.end();
				tIBFArbitraryOffset += timer.getTime();
				return false;
			}
		}
	}
	delete[] schedulable;
	timer.end();
	tIBFArbitraryOffset += timer.getTime();
	return true;
}

double SchedulabilityAnalysis::get_ibf_arbitrary_offset(Stateflow** sfs, int i, int t, int tprim) {
	double total = 0;
	double dbfi = sfs[i]->get_dbf(t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double ibfj = sfs[j]->get_ibf(tprim);
		total += ibfj;

		//if (total > tprim) return false;
	}

	return total;
}