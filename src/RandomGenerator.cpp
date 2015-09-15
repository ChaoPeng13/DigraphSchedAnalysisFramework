#include "RandomGenerator.h"
#include "GraphAlgorithms.h"

// Initialize Digraph Generating Parameters
int	RandomGenerator::DIGRAPH_SCALE = 1000;
int RandomGenerator::DIGRAPH_PERIOD[3][2] = {{2,4},{6,12},{5,10}};
int RandomGenerator::DIGRAPH_SCALE_FACTOR[5] = {1,2,4,5,10};

// Initialize Stateflow Generating Parameters
int RandomGenerator::STATEFLOW_SCALE = 1000;
int RandomGenerator::STATEFLOW_PERIOD[6] = {5, 10, 20, 25, 50, 100};


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
		digraph->calculate_gcd();
		digraph->calculate_linear_factor();

		// step 2: rescale wcets
		int run = 0;
		while(abs(digraph->linear_factor-util[i])>0.001) {
			double factor = util[i]/digraph->linear_factor;
			digraph->scale_wcet(factor);
			digraph->calculate_linear_factor();
			cout<<"Run "<<run<<":\tlinear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		}

		digraphs[i] = digraph;
	}

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

Stateflow** RandomGenerator::generate_stateflow_system(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		Stateflow* sf = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);

		//double* tran_util = Utility::uniformly_distributed(numTran,1);
		int k=0;
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;

			t->period = STATEFLOW_PERIOD[rand()%6]*STATEFLOW_SCALE;

			int wcet = (int) (t->period*(rand()%100/100.0));
			if (wcet == 0) wcet = rand()%5+1;
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
		while (abs(sf->lfac-util[i])>0.001) {
			double factor = util[i]/sf->lfac;
			sf->scale_wcet(factor);
			sf->generate_rbfs();
			sf->generate_exec_req_matrix();
			//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
			sf->calculate_linear_factor();

			//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util[i]<<"\tfact="<<sf->lfac<<endl;
		}

		sfs[i] = sf;
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


