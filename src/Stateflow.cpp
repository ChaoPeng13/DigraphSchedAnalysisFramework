#include "Stateflow.h"

#include "MaxPlusAlgebra.h"
#include "GraphAlgorithms.h"

#include "Timer.h"

double Stateflow::tDiff = 0;
double Stateflow::tCalWithoutPeriodicity = 0;
double Stateflow::tCalWithPeriodicity = 0;

Stateflow::~Stateflow() {
	//cout<<"Stateflow Destruction Start"<<endl;
	state_index.clear();
	index_state.clear();
	tran_index.clear();
	index_tran.clear();
	// release states
	for (vector<State*>::iterator iter = states.begin(); iter != states.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	states.clear();
		
	// release transitions
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	trans.clear();

	rbf_time_instances.clear();
	index_time.clear();
	time_index.clear();

	// delete rbfs
	if (rbfs != NULL) {
		for (int i=0; i<n_state; i++) {
			for (int j=0; j<n_state; j++) {
				for (int s=0; s< n_time_instance; s++) {
					delete[] rbfs[i][j][s];
				}
				delete[] rbfs[i][j];
			}
			delete[] rbfs[i];
		}
		delete[] rbfs;
	}

	//rbf_hyperperiod.clear();
	//ibf_hyperperiod.clear();
		
	// delete execution request matrix
	for(map<int, double**>::iterator iter = exec_req_matrix_power.begin(); iter != exec_req_matrix_power.end(); iter++) {
		double** temp = iter->second;
		for (int i=0; i<n_state; i++) {
			delete[] temp[i];
		}
		delete[] temp;
	}
	exec_req_matrix_power.clear();

	if (simple_digraph != NULL) delete simple_digraph;
	if (precise_digraph != NULL) delete precise_digraph;
	//if (exec_digraph != NULL) delete exec_digraph;

	// release rbf_vec
	for (vector<StartFinishTime*>::iterator iter = rbf_vec.begin(); iter != rbf_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	rbf_vec.clear();

	// release ibf_vec
	for (vector<StartFinishTime*>::iterator iter = ibf_vec.begin(); iter != ibf_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	ibf_vec.clear();

	// release abstract request function trees
	for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mmIter = mmARFT.begin(); mmIter != mmARFT.end(); mmIter++) {
		map<int,AbstractRequestFunctionTree> temp = mmIter->second;
		for (map<int,AbstractRequestFunctionTree>::iterator mIter = temp.begin(); mIter != temp.end(); mIter++) {
			AbstractRequestFunctionTree& arft = mIter->second;
			arft.root = NULL;
			while(!arft.cur_queue.empty()) arft.cur_queue.pop();
			for (vector<RFNode*>::iterator iter = arft.record.begin(); iter != arft.record.end(); iter++) {
				delete *iter;
			}
			arft.record.clear();
		}
	}

	//cout<<"Staeflow End"<<endl;
}

void Stateflow::add_state(State* state) {
	states.push_back(state);
	state_index[state] = iState;
	index_state[iState] = state;
	iState++;
}

void Stateflow::add_transition(Transition* tran) {
	trans.push_back(tran);
	tran_index[tran] = iTran;
	index_tran[iTran] = tran;

	tran->src->out.push_back(tran);
	tran->snk->in.push_back(tran);

	iTran++;
}

bool Stateflow::less_compare_trans(const Transition* t1, const Transition* t2) {
	return t1->priority < t2->priority;
}

void Stateflow::prepare() {
	// TODO:
}

void Stateflow::calculate_gcd() {
	int _gcd = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_gcd = Utility::math_gcd(_gcd, tran->period);
	}

	gcd = _gcd;
}

void Stateflow::calculate_t_gcd() {
	int _gcd = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_gcd = Utility::math_gcd(_gcd, tran->period);
		_gcd = Utility::math_gcd(_gcd, tran->wcet);
	}
	t_gcd = _gcd;
}

void Stateflow::calculate_hyperperiod() {
	int _lcm = 1;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_lcm = Utility::math_lcm(_lcm, tran->period);
	}
	hyperperiod = _lcm;
}

void Stateflow::calculate_deadlines() {
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		State* snk = (*iter)->snk;
		int deadline = (*iter)->period;
		for (list<Transition*>::iterator iter2 = snk->out.begin(); iter2 != snk->out.end(); iter2++) {
			deadline = min(deadline,(*iter2)->period);
		}
		(*iter)->deadline = deadline;
	}
}

void Stateflow::set_state_number() {
	n_state = states.size();
}

void Stateflow::generate_rbf_time_instances() {
	rbf_time_instances.insert(0);

	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		for (int i=tran->period; i<hyperperiod; i+=tran->period)
			rbf_time_instances.insert(i);
	}

	rbf_time_instances.insert(hyperperiod);
	int index = 0;
	for (set<int>::iterator iter = rbf_time_instances.begin(); iter != rbf_time_instances.end();  iter++) {
		index_time[index] = *iter;
		time_index[*iter] = index;
		index++;
	}

	n_time_instance = rbf_time_instances.size();
}

void Stateflow::generate_ibf_time_instances() {
	// TODO: It needs to think more
}

void Stateflow::generate_rbfs() {
	// Sort the transitions by the priorities
	//sort(trans.begin(),trans.end(),less_compare_trans);
	if (rbfs == NULL) { // allocate memory for rbfs
		rbfs = new double***[n_state];
		for (int i=0; i<n_state; i++) rbfs[i] = new double**[n_state];
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) rbfs[i][j] = new double*[n_time_instance];
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) for (int k=0; k<n_time_instance; k++) rbfs[i][j][k] = new double[n_time_instance];
	}

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++)
		for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
			if (i==j && s<=f) rbfs[i][j][s][f] = 0;
			else rbfs[i][j][s][f] = NEG_INFINITY;
		}
	
	// Dynamic Programming
	for (int len=1; len<n_time_instance; len++) {
		
		for (int s=0; s<n_time_instance-len; s++) {
			int f = s+len; // index of finish time
			int prev_ft = index_time[f-1]; // previous finish time

			//initialize
			for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
				if (len==1) {
					if (i==j) rbfs[i][j][s][f]=0;
				} else rbfs[i][j][s][f] = rbfs[i][j][s][f-1];
			}

			// calculate rbfs[i][j][s][f] based on the assumption we have known all the rbf[x][y][s][s..f-1].

			for (vector<Transition*>::iterator t_iter = trans.begin(); t_iter != trans.end(); t_iter++) {
				Transition* tran = *t_iter;
				int srcIndx = state_index[tran->src];
				int snkIndx = state_index[tran->snk];

				if(prev_ft%tran->period == 0) { // this transition might happen
					for (vector<State*>::iterator s_iter = states.begin(); s_iter != states.end(); s_iter++) {
						State* state = *s_iter;
						int sIndx = state_index[state];

						if (len == 1 && sIndx == srcIndx)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], (double)tran->wcet);
						else if (len >1)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], rbfs[sIndx][srcIndx][s][f-1]+tran->wcet);
					}
				}
			}
		}
	}

}

void Stateflow::generate_ibfs() {
	// TODO: It needs to think more
}

void Stateflow::generate_exec_req_matrix() {
	if (exec_req_matrix == NULL) { // allocate memory
		exec_req_matrix = new double*[n_state];
		for (int i=0; i<n_state; i++) exec_req_matrix[i] = new double[n_state];
	}

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		exec_req_matrix[i][j] = rbfs[i][j][0][n_time_instance-1];
	}
	exec_req_matrix_power[1] = exec_req_matrix;
}

void Stateflow::calculate_exec_req_matrix_power(int tf) {
	Timer timer;
	/// Calculate all the execution request matrix from 2 to tf
	for (int t=2; t<tf; t++) {
		if (exec_req_matrix_power.find(t) != exec_req_matrix_power.end())
			continue;
		double** temp;
		if (t>gdef+gper && isIrred && exec_req_matrix_power.find(t-gper) != exec_req_matrix_power.end()) {
			//cout<<"a-ha!"<<endl;
			timer.start();
			temp = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(exec_req_matrix_power[t-gper],n_state,lfac*hyperperiod,gper);
			timer.end();
			double t0 = timer.getTime();

			timer.start();
			//double** temp1 = MaxPlusAlgebra::multiply_maxplus_matrix(exec_req_matrix_power[t-1],exec_req_matrix,n_state);
			double** temp1 = MaxPlusAlgebra::power_maxplus_matrix(exec_req_matrix,n_state, t);
			timer.end();
			double t1 = timer.getTime();

			for (int i=0; i<n_state; i++) delete[] temp1[i];
			delete[] temp1;

			cout<<t0<<"\t"<<t1<<"\t"<<t1-t0<<endl;

			tDiff += t1-t0;
			tCalWithoutPeriodicity += t1;
			tCalWithPeriodicity += t0;
		}
		else
			temp = MaxPlusAlgebra::multiply_maxplus_matrix(exec_req_matrix_power[t-1],exec_req_matrix,n_state);
		exec_req_matrix_power[t] = temp;
	}
}

void Stateflow::calculate_generialized_factor() {
	gfac = Utility::creat_matrix(n_state,n_state);
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) gfac[i][j] = lfac;
}

void Stateflow::calculate_generialized_period() {
	if (!isIrred) return;
	gper = MaxPlusAlgebra::calculate_linear_period(exec_req_matrix,n_state,lfac*hyperperiod,false);
}

void Stateflow::calculate_generialized_defect() {
	if (!isIrred) 
		gdef = 100*scale;
	else
		gdef = MaxPlusAlgebra::calculate_linear_defect(exec_req_matrix_power, n_state, lfac*hyperperiod, gper, tf0);
}

void Stateflow::check_irreducible() {
	//this->lfac = precise_digraph->linear_factor;
	this->isIrred = simple_digraph->strongly_connected;
}

void Stateflow::calculate_linear_factor() {
	if (isIrred)
		lfac = GraphAlgorithms::calculate_maximum_cycle_mean(exec_req_matrix, n_state)/hyperperiod;
	else {
		exec_digraph = new GeneralDirectedGraph(n_state, exec_req_matrix);
		//write_graphviz(cout);
		//exec_digraph->write_graphviz(cout);
		exec_digraph->generate_strongly_connected_components();
		exec_digraph->calculate_untilization();
		lfac = exec_digraph->util/hyperperiod;
		delete exec_digraph;
	}
}

void Stateflow::scale_wcet(double factor) {
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* t = *iter;
		int wcet = (int)(factor*t->wcet);
		if (wcet == 0) wcet = rand()%5+1; // force the wcet not to be 0
		t->wcet = wcet;
	}
}

void Stateflow::calculate_csum() {
	int _csum = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		//int temp = tran->wcet*(hyperperiod/tran->period);
		_csum += tran->wcet*(hyperperiod/tran->period);
	}
	csum = _csum;
}

void Stateflow::calculate_tf0(Stateflow** stateflows, int i) {
	double util = 0;
	double sum = 0;
	for (int k=0; k<=i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->csum*2;
		if (k==i) sum -= sf->csum;
	}

	if (util >= 1) 
		tf0 = POS_INFINITY;
	else 
		tf0 = sum/(1-util);
}

/// Generate a simple digraph for stateflow
void Stateflow::generate_simple_digraph() {
	if (simple_digraph != NULL) return;

	simple_digraph = GraphAlgorithms::generate_simple_digraph(this);
	// prepare for calculating linear upper bounds
	simple_digraph->prepare_digraph();
}

/// Generate a precise digraph for stateflow
void Stateflow::generate_precise_digraph() {
	if (precise_digraph != NULL) return;
	precise_digraph = GraphAlgorithms::generate_precise_digraph(this);

	// prepare for calculating linear upper bounds
	precise_digraph->prepare_digraph();

	//lfac = precise_digraph->linear_factor;
}

/// Calculate the tigher linear upper bounds
/// true: operating on precise digraph
/// false: operating on simple digraph
void Stateflow::calculate_linear_upper_bounds(bool precise) {
	if(precise) {
		// generate_precise_digraph();
		crbf = precise_digraph->c_rbf;
		cibf = precise_digraph->c_ibf;
		cdbf = precise_digraph->c_dbf;
	} else {
		//  generate_simple_digraph();
		crbf = simple_digraph->c_rbf;
		cdbf = simple_digraph->c_dbf;
		cibf = simple_digraph->c_ibf;
	}
}

void Stateflow::calculate_tf1(Stateflow** stateflows, int i) {
	double util = stateflows[i]->lfac;
	double sum = stateflows[i]->cdbf;
	for (int k=0; k<i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->crbf;
	}

	if (util >= 1) 
		tf1 = POS_INFINITY;
	else
		tf1 = sum/(1-util);
}

void Stateflow::calculate_tf2(Stateflow** stateflows, int i) {
	double util = stateflows[i]->lfac;
	double sum = stateflows[i]->cdbf;
	for (int k=0; k<i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->cibf;
	}

	if (util >= 1) 
		tf2 = POS_INFINITY;
	else
		tf2 = sum/(1-util);
}

/// Static offset
double Stateflow::get_rbf(int start, int finish) { // return rbf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int fprim = *rbf_time_instances.lower_bound(f);
	int nf = finish/hyperperiod;
	int n = nf - ns;
	f = n*hyperperiod+fprim;

	for (vector<StartFinishTime*>::iterator iter = rbf_vec.begin(); iter != rbf_vec.end(); iter++) {
		StartFinishTime* sft = *iter;
		if (sft->start == s && sft->finish == f)
			return sft->value;
	}
	
	double rt = 0;

	if (f<=hyperperiod) // s and f are in the same hyperperiod
		rt = calculate_rbf_within_one_hyperperiod(s,f);
	else 
		rt = calculate_rbf_within_multiple_hyperperiods(s,f);

	StartFinishTime* sft = new StartFinishTime(s,f,rt);
	rbf_vec.push_back(sft);

	return rt;
}

double Stateflow::calculate_rbf_within_one_hyperperiod(int start, int finish) {
	double rt = NEG_INFINITY;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		rt = max(rt, calculate_rbf_within_one_hyperperiod(i,j,start,finish));
	}
	return rt;
}

double Stateflow::calculate_rbf_within_one_hyperperiod(int i, int j, int start, int finish) {
	int sIndx = time_index[start];
	int fIndx = time_index[finish];

	return rbfs[i][j][sIndx][fIndx];
}

double Stateflow::calculate_rbf_within_multiple_hyperperiods(int start, int finish) {

	int n = finish/hyperperiod;
	int f = finish - n*hyperperiod;

	double** prev = Utility::creat_matrix(n_state,n_state);
	double** midd = Utility::creat_matrix(n_state,n_state);
	double** post = Utility::creat_matrix(n_state,n_state);

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		prev[i][j] = 0;
		midd[i][j] = 0;
		post[i][j] = 0;
	}
	

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		prev[i][j] = calculate_rbf_within_one_hyperperiod(i,j,start,hyperperiod);
	}
	

	if (n>1) {
		if (exec_req_matrix_power.find(n)==exec_req_matrix_power.end()) 
			calculate_exec_req_matrix_power(n);
		
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state;j++)
			midd[i][j] = exec_req_matrix_power[n-1][i][j];
	}
	
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		post[i][j] = calculate_rbf_within_one_hyperperiod(i,j,0,f);
	}

	double rbf = 0;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		double rbfij = 0;
		if (n!=1) {
			for (int k=0; k<n_state; k++) for (int l=0; l<n_state; l++) {
				rbfij = max(rbfij, prev[i][k]+midd[k][l]+post[l][j]);
			}
		}
		else {// n==1, we should not use the middle matrix in case of the errors of states
			for (int k=0; k<n_state; k++) {
				rbfij = max(rbfij, prev[i][k]+post[k][j]);
			}
		}
		rbf = max(rbf, rbfij);
	}

	// release prev, midd and post matrices
	for (int i=0; i<n_state; i++) {
		delete[] prev[i];
		delete[] midd[i];
		delete[] post[i];
	}
	delete[] prev;
	delete[] midd;
	delete[] post;

	return rbf;
}

double Stateflow::get_ibf(int start, int finish) { // return ibf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int nf = finish/hyperperiod;
	int n = nf - ns;
	f = n*hyperperiod+f;

	for (vector<StartFinishTime*>::iterator iter = ibf_vec.begin(); iter != ibf_vec.end(); iter++) {
		StartFinishTime* sft = *iter;
		if (sft->start == s && sft->finish == f)
			return sft->value;
	}
	
	double rt = 0;

	if (f<=hyperperiod) // s and f are in the same hyperperiod
		rt = calculate_ibf_within_one_hyperperiod(s,f);
	else 
		rt = calculate_ibf_within_multiple_hyperperiods(s,f);

	StartFinishTime* sft = new StartFinishTime(s,f,rt);
	ibf_vec.push_back(sft);

	return rt;
}

double Stateflow::calculate_ibf_within_one_hyperperiod(int start, int finish) {
	double rt = NEG_INFINITY;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		rt = max(rt, calculate_ibf_within_one_hyperperiod(i,j,start,finish));
	}
	return rt;
}

double Stateflow::calculate_ibf_within_one_hyperperiod(int i, int j, int start, int finish) {
	int sIndx = time_index[start];
	int fIndx = 0;

	int fh = *rbf_time_instances.lower_bound(finish);
	if (fh == finish) {
		fIndx = time_index[finish];
		return rbfs[i][j][sIndx][fIndx];
	}
	
	int fl = *(--rbf_time_instances.lower_bound(finish));
	int fIndx0 = time_index[fl];

	double ibfij = rbfs[i][j][sIndx][fIndx0];
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		int src = state_index[tran->src];
		int snk = state_index[tran->snk];
		if (snk != j) continue;
		if (fl%tran->period==0) {
			double dmin = min(tran->wcet, finish-fl);
			double rbfij = rbfs[i][src][sIndx][fIndx0];
			if (rbfij == NEG_INFINITY) continue;
			ibfij = max(ibfij, rbfij+dmin);
		}
	}

	return ibfij;
}

double Stateflow::calculate_ibf_within_multiple_hyperperiods(int start, int finish) {

	int n = finish/hyperperiod;
	int f = finish - n*hyperperiod;

	double** prev = Utility::creat_matrix(n_state,n_state);
	double** midd = Utility::creat_matrix(n_state,n_state);
	double** post = Utility::creat_matrix(n_state,n_state);

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		prev[i][j] = 0;
		midd[i][j] = 0;
		post[i][j] = 0;
	}
	

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		prev[i][j] = calculate_rbf_within_one_hyperperiod(i,j,start,hyperperiod);
	}
	

	if (n>1) {
		if (exec_req_matrix_power.find(n)==exec_req_matrix_power.end()) 
			calculate_exec_req_matrix_power(n);

		for (int i=0; i<n_state; i++) for (int j=0; j<n_state;j++)
			midd[i][j] = exec_req_matrix_power[n-1][i][j];
	}

	
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		post[i][j] = calculate_ibf_within_one_hyperperiod(i,j,0,f);
	}

	double ibf = 0;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		double ibfij = 0;
		if (n!=1) {
			for (int k=0; k<n_state; k++) for (int l=0; l<n_state; l++) {
				ibfij = max(ibfij, prev[i][k]+midd[k][l]+post[l][j]);
			}
		}
		else {// n==1, we should not use the middle matrix in case of the errors of states
			for (int k=0; k<n_state; k++) {
				ibfij = max(ibfij, prev[i][k]+post[k][j]);
			}
		}
		ibf = max(ibf, ibfij);
	}

	// release prev, midd and post matrices
	for (int i=0; i<n_state; i++) {
		delete[] prev[i];
		delete[] midd[i];
		delete[] post[i];
	}
	delete[] prev;
	delete[] midd;
	delete[] post;

	return ibf;
}

/// Note: dbf[s,f] = rbf[s,lower(f))
double Stateflow::get_dbf(int start, int finish) { // return dbf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int hf = *rbf_time_instances.lower_bound(f);
	int hl = 0;
	if (f!=0) hl = *(--rbf_time_instances.lower_bound(f));
	int nf = finish/hyperperiod;

	int n = nf - ns;

	if (f == hf) {
		f = n*hyperperiod+hf;
	} else if (f-hl >= gcd) {
		f = n*hyperperiod+hf;
	} else {
		f = n*hyperperiod+hl;
	}
	
	if (s > f) return NEG_INFINITY;
	if (s == f) return 0;

	return get_rbf(s,f);

}

/// Arbitrary offset
double Stateflow::get_rbf(int t) { // return rbf(t)
	if (t==0) return 0;

	double rbf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double rbfsf = get_rbf(s,f);
		rbf = max(rbf, rbfsf);
	}
	return rbf;
}

double Stateflow::get_ibf(int t) { // return ibf(t)
	if (t==0) return 0;
	
	double ibf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double ibfsf = get_ibf(s,f);
		ibf = max(ibf, ibfsf);
	}
	return ibf;
}

double Stateflow::get_dbf(int t) { // return dbf(t)
	if (t==0) return 0;

	double dbf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double dbfsf = get_dbf(s,f);
		dbf = max(dbf, dbfsf);
	}
	return dbf;
}

void Stateflow::write_graphviz(ostream& out) {
	out << "digraph G {" <<endl;
	for (vector<State*>::iterator iter = states.begin(); iter != states.end(); iter++) {
		State* state = *iter;
		out << state->name << " [label=\" " << state->name << " \"]" << endl;
	}
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* transition = *iter;
		out << transition->src->name << " -> " << transition->snk->name 
			<< " [label=\" " << transition->wcet << " / " << transition->period << " \"]" << endl;
	}
	out << "}" <<endl;
}


