#include "RandomGenerator.h"
#include "GraphAlgorithms.h"
#include "FileReader.h"

// Initialize Digraph Generating Parameters
int	RandomGenerator::DIGRAPH_SCALE = 1000;
int RandomGenerator::DIGRAPH_PERIOD[3][2] = {{2,4},{6,12},{5,10}};
int RandomGenerator::DIGRAPH_SCALE_FACTOR[5] = {1,2,4,5,10};

// Initialize Stateflow Generating Parameters
int RandomGenerator::STATEFLOW_SCALE = 1000;
int RandomGenerator::STATEFLOW_PERIOD[6] = {5, 10, 20, 25, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD2[5] = {10, 20, 25, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD3[4] = {10, 20, 25, 50};
int RandomGenerator::STATEFLOW_PERIOD4[5] = {5, 10, 20, 25, 50};
//int RandomGenerator::STATEFLOW_FACTOR[6] = {1, 2, 4, 5, 10, 20};
int RandomGenerator::STATEFLOW_BASE[8] = {1100, 1300,2300, 3100, 3700, 3000, 5000, 7000};
//int RandomGenerator::STATEFLOW_BASE[5] = {10000, 11000, 13000, 15000, 17000};
//int RandomGenerator::STATEFLOW_BASE[8] = {110, 130,230, 310, 370, 300, 500, 700};
int RandomGenerator::STATEFLOW_FACTOR[5] = {1, 2, 4, 5, 10};
int RandomGenerator::STATEFLOW_FACTOR2[4] = {1, 2, 4, 5};
int RandomGenerator::STATEFLOW_FACTOR3[3] = {1, 2, 4};
int RandomGenerator::STATEFLOW_FACTOR4[2] = {1, 2};

Digraph* RandomGenerator::generate_one_digraph(int index, int scale, int nNode, int nEdge) {
	Digraph* digraph = new Digraph(index, scale);

	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	connected.push_back(nodes[rand()%nNode]);
	for (int i=0; i<nNode; i++) {
		Node* newNode = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == newNode) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		Node* connectedNode = connected.at(rand()%connected.size());
		Edge* edge = NULL;

		int new2Connected = rand()%2;

		if (new2Connected) 
			edge = new Edge(newNode, connectedNode);
		else 
			edge = new Edge(connectedNode, newNode);

		digraph->add_edge(edge);
		connected.push_back(newNode);
	}

	for (int i= nNode-1; i<nEdge; i++) {
		Node* src = nodes[rand()%nNode];
		Node* snk = nodes[rand()%nNode];

		Edge* edge = new Edge(src,snk);
		digraph->add_edge(edge);
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		for (int i=0; i<nEdge; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Edge* edge = digraph->edge_vec.at(i);

				// reverse the direction of the edge
				edge->src_node->in.push_back(edge);
				edge->src_node->out.remove(edge);

				edge->snk_node->out.push_back(edge);
				edge->snk_node->in.remove(edge);

				Node* temp = edge->src_node;
				edge->src_node = edge->snk_node;
				edge->snk_node = temp;
			}
		}
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph(int index, int scale, int nNode) {
	Digraph* digraph = new Digraph(index, scale);

	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	return digraph;
}

Digraph** RandomGenerator::generate_digraph_system(int num, int maxNode, double tUtil) {
	Digraph** digraphs = new Digraph*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numNode = rand()%maxNode+1;

		Digraph* digraph = generate_one_digraph(i,DIGRAPH_SCALE,numNode);

		int base = calculate_base();

		// generate random edges with random separation time
		for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
			Edge* edge = *iter;

			int separationTime = calculate_separation_time(base);
			edge->set_separation_time(separationTime*DIGRAPH_SCALE);
		}

		// generate random nodes with random deadline and wcet
		for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
			Node* node = *iter;

			int minSeparationTime = INT_MAX;
			for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) 
				minSeparationTime = min(minSeparationTime, (*e)->separationTime);

			int deadline = minSeparationTime;
			node->set_deadline(deadline);

			int wcet = rand()%(deadline/2)+1;
			//if (wcet == 0) wcet = rand()%5+1;
			node->set_wcet(wcet);
		}
		
		/// scale wcets to the expected utilization
		// step 1: calculate the linear factor
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_period_gcd();
		digraph->calculate_linear_factor();
		//cout<<"linear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		// step 2: rescale wcets
		int run = 0;
		while(abs(digraph->linear_factor-util[i])>0.001) {
			double factor = util[i]/digraph->linear_factor;
			digraph->scale_wcet(factor);
			digraph->calculate_linear_factor();
			//cout<<"Run "<<run<<":\tlinear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		}

		digraphs[i] = digraph;
	}

	// reorder digraphs according to the utilization (non-decreasing order)
	// It is possible to considered more.
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		order[i] = 0;
		for (int j=0; j<num; j++) {
			Digraph* hdg = digraphs[j];
			if (hdg->linear_factor < digraph->linear_factor) order[i]++;
		}
	}

	Digraph** copy = new Digraph*[num];
	for (int i=0; i<num; i++) copy[i] = digraphs[i];
	for (int i=0; i<num; i++) digraphs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return digraphs;
}


int RandomGenerator::calculate_base() {
	int base = 1;
	int num_ext = rand()%3+1;
	for (int j=0; j<num_ext; j++) base *= DIGRAPH_PERIOD[rand()%3][rand()%2];

	return base;
}

int RandomGenerator::calculate_separation_time(int base) {
	int factor = DIGRAPH_SCALE_FACTOR[rand()%5];

	return factor*base;
}

Stateflow* RandomGenerator::generate_one_stateflow(int index, int scale, int nState, int nTran) {
	Stateflow* sf = new Stateflow(index, scale);

	State** states = new State*[nState];

	for (int i=0; i<nState; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		State* s = new State(name, i);
		states[i] = s;
		sf->add_state(s);
	}

	vector<State*> connected;
	connected.push_back(states[rand()%nState]);
	for (int i=0; i<nState; i++) {
		State* s = states[i];
		
		bool isContained = false;
		for (vector<State*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == s) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		State* s2 = connected.at(rand()%connected.size());
		Transition* t = NULL;

		int b = rand()%2;

		if (b)
			t = new Transition(s,s2);
		else 
			t = new Transition(s2,s);

		sf->add_transition(t);
		connected.push_back(s);
	}

	for (int i=nState-1; i<nTran; i++) {
		State* src = states[rand()%nState];
		State* snk = states[rand()%nState];

		Transition* t = new Transition(src,snk);
		sf->add_transition(t);
	}
	//sf->set_state_number();
	//sf->write_graphviz(cout);
	//int count = 0;
	while (!GraphAlgorithms::isCyclic(sf)) {
		for (int i=0; i<nTran; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Transition* t = sf->trans.at(i);
				
				// reverse the direction of the transition
				t->src->in.push_back(t);
				t->src->out.remove(t);

				t->snk->out.push_back(t);
				t->snk->in.remove(t);

				State* temp = t->src;
				t->src = t->snk;
				t->snk = temp;
			}
		}
		//count ++;
	}
	//cout<<count<<endl;
	//sf->write_graphviz(cout);
	connected.clear();
	delete[] states;
	return sf;
}

Stateflow* RandomGenerator::generate_one_stateflow2(int index, int scale, int nState, int nTran) {
	Stateflow* sf = new Stateflow(index, scale);

	State** states = new State*[nState];

	for (int i=0; i<nState; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		State* s = new State(name, i);
		states[i] = s;
		sf->add_state(s);
	}

	vector<State*> connected;
	connected.push_back(states[rand()%nState]);
	for (int i=0; i<nState; i++) {
		State* s = states[i];
		
		bool isContained = false;
		for (vector<State*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == s) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		State* s2 = connected.at(rand()%connected.size());
		Transition* t = NULL;

		int b = rand()%2;

		if (b)
			t = new Transition(s,s2);
		else 
			t = new Transition(s2,s);

		sf->add_transition(t);
		connected.push_back(s);
	}

	for (int i=nState-1; i<nTran; i++) {
		State* src = states[rand()%nState];
		int count = 0;
		bool exceed = false;
		while (src->out.size() >= 2) {
			src = states[rand()%nState];
			count ++;
			if (count > 100) {
				exceed = true;
				break;
			}
		}
		if (exceed) {nTran = i+1; break;}

		State* snk = states[rand()%nState];

		Transition* t = new Transition(src,snk);
		sf->add_transition(t);
	}
	//sf->set_state_number();
	//sf->write_graphviz(cout);
	//int count = 0;
	while (!GraphAlgorithms::isCyclic(sf)) {
		for (int i=0; i<nTran; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Transition* t = sf->trans.at(i);
				if (t->snk->out.size() >= 1) continue;
				
				// reverse the direction of the transition
				t->src->in.push_back(t);
				t->src->out.remove(t);

				t->snk->out.push_back(t);
				t->snk->in.remove(t);

				State* temp = t->src;
				t->src = t->snk;
				t->snk = temp;
			}
		}
		//count ++;
	}
	//cout<<count<<endl;
	//sf->write_graphviz(cout);
	connected.clear();
	delete[] states;
	return sf;
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_0(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
	int period_base = STATEFLOW_PERIOD3[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR4)/sizeof(STATEFLOW_FACTOR4[0]);
		int period_factor = STATEFLOW_FACTOR4[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_1(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD4)/sizeof(STATEFLOW_PERIOD4[0]);
	int period_base = STATEFLOW_PERIOD4[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_2(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD4)/sizeof(STATEFLOW_PERIOD4[0]);
	int period_base = STATEFLOW_PERIOD4[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR4)/sizeof(STATEFLOW_FACTOR4[0]);
		int period_factor = STATEFLOW_FACTOR4[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_3(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
	int period_base = STATEFLOW_PERIOD3[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_0(Stateflow* sf, int scale) {
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
		t->period = STATEFLOW_PERIOD[rand()%length]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_1(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
	int period_base = STATEFLOW_PERIOD[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR2)/sizeof(STATEFLOW_FACTOR2[0]);
		int period_factor = STATEFLOW_FACTOR2[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_2(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
	int period_base = STATEFLOW_PERIOD[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_wcets_for_one_stateflow(Stateflow* sf, double util) {
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int wcet = (int) (t->period*(rand()%50/50.0+0.5));
		if (wcet == 0) wcet = rand()%5+1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.0001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}
}

Stateflow* RandomGenerator::generate_one_stateflow_with_util(int index, int scale, double util, int nState, int nTran) {
	Stateflow* sf = generate_one_stateflow(index,scale,nState, nTran);

	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
		t->period = STATEFLOW_PERIOD[rand()%length]*scale;

		int wcet = (int) (t->period*(rand()%100/100.0));
		if (wcet == 0) wcet = rand()%5+1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;

		t->priority = k++;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}

	return sf;
}

Stateflow* RandomGenerator::generate_one_stateflow_with_util2(int index, int scale, double util, int nState, int nTran) {
	Stateflow* sf = generate_one_stateflow(index,scale,nState, nTran);

	int base_length = sizeof(STATEFLOW_BASE)/sizeof(STATEFLOW_BASE[0]);
	int base = STATEFLOW_BASE[rand()%base_length];
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		int factor_length = sizeof(STATEFLOW_FACTOR)/sizeof(STATEFLOW_FACTOR[0]);
		int factor = STATEFLOW_FACTOR[rand()%factor_length];
		t->period = base*factor;

		//int wcet = (int) (t->period*(rand()%100/100.0));
		//if (wcet == 0) wcet = rand()%5+1;
		int wcet = factor*(rand()%5+1);
		t->wcet = wcet;

		t->priority = k++;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}

	return sf;
}

Stateflow** RandomGenerator::generate_stateflow_system(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system2(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util2(i,STATEFLOW_SCALE, util[i], numState,numTran);;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system3(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	int bound = 2*num/3;

	for (int i=0; i<bound; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util(i,STATEFLOW_SCALE, util[i], numState,numTran);
	}

	for (int i=bound; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util2(i,STATEFLOW_SCALE, util[i], numState,numTran);
	}

	delete[] util;
	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system4(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num-1; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num-1];
	for (int i=0; i<num-1; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num-1; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num-1];
	for (int i=0; i<num-1; i++) copy[i] = sfs[i];
	for (int i=0; i<num-1; i++) sfs[order[i]] = copy[i];

	double* util = Utility::uniformly_distributed(num-1, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num-1; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	string name = "Input\\OneSpecial.dot";
	const char *p = name.c_str();

	Stateflow* sf = FileReader::ReadOneStateflow(1,p);
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	sfs[num-1] = sf;

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system5(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system6(int num, int maxState, double tUtil, int period_choice) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system7(int num, int maxState, double tUtil, int period_choice) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system8(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			/*
			if (digraph->strongly_connected && r<scc_probability) {
				delete digraph;
				break;
			}
			else if (!digraph->strongly_connected && r<1-scc_probability) {
				delete digraph;
				break;
			}
			else {
				delete digraph;
				delete sf;
			}
			*/
			delete digraph;
			delete sf;

		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_exact_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		// number of states = {20%--1, 15%--2, 
		// 15%--3, 10%--4, 10%--5, 10%--6, 5%--7, 
		// 5%--8, 5%--9, 5%--10}, which makes the smaller number 
		// of states has higher proportion
		/*
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.5) numState = 3;
		else if (r<0.6) numState = 4;
		else if (r<0.7) numState = 5;
		else if (r<0.8) numState = 6;
		else if (r<0.85) numState = 7;
		else if (r<0.9) numState = 8;
		else if (r<0.95) numState = 9;
		else if (r<=1) numState = 10;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		// number of states = {30%--1, 20%--2, 20%--3, 20%--4, 10%--5}
		if (r<0.3) numState = 1;
		else if (r<0.5) numState = 2;
		else if (r<0.7) numState = 3;
		else if (r<0.9) numState = 4;
		else if (r<=1) numState = 5;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		
		int numTran = calculate_num_edge(numState);
		numTran = (numTran+numState)/2;
		cout<<"numState="<<numState<<"\t"<<"numTran"<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;

		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		/*
		// fist configure
		// number of states = {20%--1, 15%--2, 
		// 15%--3, 10%--4, 10%--5, 10%--6, 5%--7, 
		// 5%--8, 5%--9, 5%--10}, which makes the smaller number 
		// of states has higher proportion
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.5) numState = 3;
		else if (r<0.6) numState = 4;
		else if (r<0.7) numState = 5;
		else if (r<0.8) numState = 6;
		else if (r<0.85) numState = 7;
		else if (r<0.9) numState = 8;
		else if (r<0.95) numState = 9;
		else if (r<=1) numState = 10;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		
		// second configure
		// number of states = {15%--1, 15%--2, 
		// 10%--3, 10%--4, 10%--5, 5%---6, 5%---7, 
		// 5%---8, 5%---9, 5%---10, 5%---11, 4%---12, 
		// 3%---13, 2%---14, 1%---15}

		if (r<0.15) numState = 1;
		else if (r<0.30) numState = 2;
		else if (r<0.40) numState = 3;
		else if (r<0.50) numState = 4;
		else if (r<0.60) numState = 5;
		else if (r<0.65) numState = 6;
		else if (r<0.70) numState = 7;
		else if (r<0.75) numState = 8;
		else if (r<0.80) numState = 9;
		else if (r<0.85) numState = 10;
		else if (r<0.9) numState = 11;
		else if (r<0.94) numState = 12;
		else if (r<0.97) numState = 13;
		else if (r<=0.99) numState = 14;
		else if (r<=1) numState = 15;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		
		
		// third configure
		//numState = rand()%maxState+1;

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran"<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

int RandomGenerator::calculate_num_edge(int numNode) {
	int divisor = 1;
    for(int j=0; j<numNode; j++) divisor = divisor + rand()%2;
    	
    int s = ceil(1.0*(rand()%numNode)/divisor); 
    s = min(s, numNode + 1);
    	
    int numEdge = numNode;
    for(int j=0; j<s; j++) numEdge = numEdge + rand()%(numNode + 1);
		
    return numEdge;
}