#include "GraphAlgorithms.h"

#include "MaxPlusAlgebra.h"
#include <algorithm>

/// /brief an implementation of the Karp's algorithm
/// R.M. Karp, A characterization of the minimum cycle mean in a digraph, Discrete Math. 23 (1978) 309�11.
/// Note: this function is just right when operating on UDRT in case of the length.
double GraphAlgorithms::calculate_maximum_cycle_mean(Digraph* digraph) {
	digraph->calculate_period_gcd();
	UnitDigraph* unitD = new UnitDigraph(digraph);
	digraph->unit_digraph = unitD;

	if(digraph->edge_vec.size() == 0)
	{
		digraph->linear_factor = 0;
		return 0;
	}

	//unitD->calculate_gcd();
	unitD->scc_generate_unit_digraph();
	//unitD->write_graphviz(cout);

	int n = unitD->node_vec.size();
	map<Node*, int> index_map;
	int index = 0;
	for (vector<Node*>::iterator iter = unitD->node_vec.begin(); iter != unitD->node_vec.end(); iter++) {
		index_map.insert(pair<Node*,int>(*iter,index++));
	}

	// first index: length; second index: vertices;
	double** D = new double*[n+1];
	for (int i=0; i<=n; i++) D[i] = new double[n];

	for (int k=0; k<=n; k++) for (int v=0; v<n; v++) D[k][v] = NEG_INFINITY;
	for (int v=0; v<n; v++) D[0][v] = 0;

	for (int k=1; k<=n; k++) {
		for (vector<Node*>::iterator iter=unitD->node_vec.begin(); iter != unitD->node_vec.end(); iter++) {
			int v = unitD->node_to_index[*iter];
			list<Edge*> in = (*iter)->in;
			for (list<Edge*>::iterator e_iter = in.begin(); e_iter != in.end(); e_iter++) {
				Node* src = (*e_iter)->src_node;
				int srcIndx = unitD->node_to_index[src];
				D[k][v] = max(D[k][v], D[k-1][srcIndx]+src->wcet);
			}
		}
	}

	// Utility::output_matrix(D,n+1,n);

	double lambda = NEG_INFINITY;
	double* M = new double[n];
	for (int v=0; v<n; v++) {
		M[v] = POS_INFINITY;
		for (int k=0; k<n; k++) M[v] = min(M[v], (D[n][v]-D[k][v])/(n-k));

		lambda = max(lambda, M[v]);
	}

	// release
	for (int i=0; i<=n; i++) delete[] D[i];
	delete[] D;

	delete[] M;


	lambda = lambda/digraph->pGCD;
	digraph->linear_factor = lambda;
	return lambda;
}

/*
 * \brief An implementation of the iterative algorithm for calculating utilization
 * in Stigge et al. The digraph real-time task model. RTAS2011.
 * In case of the large number of nodes in the unit digraph.
 */
double GraphAlgorithms::calculate_maximum_cycle_mean2(Digraph* digraph) {
	vector<vector<UtilizationTriple*>> UT;
	// initial UT0
	vector<UtilizationTriple*> UT0;
	for ( int i=0; i<digraph->node_vec.size(); i++) {
		Node* node = digraph->node_vec.at(i);
		UtilizationTriple* ut = new UtilizationTriple(0,0,node,node);
		UT0.push_back(ut);
	}
	UT.push_back(UT0);

	for (int i=0; i<digraph->node_vec.size(); i++) {
		vector<UtilizationTriple*> UTk;
		vector<UtilizationTriple*> PrevUT = UT.back();

		for (vector<UtilizationTriple*>::iterator iter = PrevUT.begin(); iter != PrevUT.end(); iter++) {
			UtilizationTriple* ut = *iter;
			Node* start = ut->start;
			Node* end = ut->end;
			for (list<Edge*>::iterator eIter = end->out.begin(); eIter != end->out.end(); eIter++) {
				Edge* edge = *eIter;
				Node* snk = edge->snk_node;

				int e = ut->e + end->wcet;
				int p = ut->p + edge->separationTime;

				bool dominated = false;
				for (vector<UtilizationTriple*>::iterator it = UTk.begin(); it != UTk.end(); it++) {
					UtilizationTriple* ut2 = *it;
					if (ut2->e >= e && ut2->p <= p && start == ut2->start  && snk == ut2->end) {
						dominated = true;
						break;
					}
					if (ut2->e <= e && ut2->p >= p && start == ut2->start  && snk == ut2->end) {
						ut2->e = e;
						ut2->p = p;
						dominated = true;
						break;
					}
				}

				if (!dominated) {
					UtilizationTriple* nut = new UtilizationTriple(e,p,start,snk);
					UTk.push_back(nut);
				}
			}
		}
		if (!UTk.empty())
			UT.push_back(UTk);
	}

	double max = 0;
	for (vector<vector<UtilizationTriple*>>::iterator iter = UT.begin(); iter != UT.end(); iter++) {
		vector<UtilizationTriple*> UTk = *iter;
		for (vector<UtilizationTriple*>::iterator it = UTk.begin(); it != UTk.end(); it++) {
			UtilizationTriple* ut = *it;
			if (ut->start == ut->end)
				max = std::max(max,1.0*ut->e/ut->p);
		}
	}

	return max;
}

/// /brief an implementation of the Karp's algorithm
/// R.M. Karp, A characterization of the minimum cycle mean in a digraph, Discrete Math. 23 (1978) 309�11.
/// Note: this function is just right when operating on UDRT in case of the length.
double GraphAlgorithms::calculate_maximum_cycle_mean(GeneralDirectedGraph* gDigraph) {
	if(gDigraph->edge_vec.size() == 0)
	{
		gDigraph->util = 0;
		return 0;
	}

	int n = gDigraph->node_vec.size();
	map<Node*, int> node_to_index;
	int index = 0;
	for (vector<Node*>::iterator iter = gDigraph->node_vec.begin(); iter != gDigraph->node_vec.end(); iter++) {
		node_to_index[*iter]=index++;
	}

	// first index: length; second index: vertices;
	double** D = new double*[n+1];
	for (int i=0; i<=n; i++) D[i] = new double[n];

	for (int k=0; k<=n; k++) for (int v=0; v<n; v++) D[k][v] = NEG_INFINITY;
	for (int v=0; v<n; v++) D[0][v] = 0;

	for (int k=1; k<=n; k++) {
		for (vector<Node*>::iterator iter= gDigraph->node_vec.begin(); iter != gDigraph->node_vec.end(); iter++) {
			int v = node_to_index[*iter];
			list<Edge*> in = (*iter)->scc_in;
			for (list<Edge*>::iterator e_iter = in.begin(); e_iter != in.end(); e_iter++) {
				Node* src = (*e_iter)->src_node;
				int srcIndx = node_to_index[src];
				D[k][v] = max(D[k][v], D[k-1][srcIndx]+(*e_iter)->weight);
			}
		}
	}

	// Utility::output_matrix(D,n+1,n);

	double lambda = NEG_INFINITY;
	double* M = new double[n];
	for (int v=0; v<n; v++) {
		M[v] = POS_INFINITY;
		for (int k=0; k<n; k++) M[v] = min(M[v], (D[n][v]-D[k][v])/(n-k));

		lambda = max(lambda, M[v]);
	}
	gDigraph->util = lambda;

	// release D
	for (int i=0; i<=n; i++)
		delete[] D[i];
	delete[] D;

	// release M
	delete[] M;

	return lambda;
}

double GraphAlgorithms::calculate_maximum_cycle_mean(double** A, int n) {
	// first index: length; second index: vertices;
	double** D = new double*[n+1];
	for (int i=0; i<=n; i++) D[i] = new double[n];

	for (int k=0; k<=n; k++) for (int v=0; v<n; v++) D[k][v] = NEG_INFINITY;
	D[0][0] = 0;

	for (int k=1; k<=n; k++) {
		for (int v=0; v<n; v++) for (int v1=0; v1<n; v1++) {
			if (A[v1][v] > 0) D[k][v] = max(D[k][v], D[k-1][v1]+A[v1][v]);
		}
	}

	double lambda = NEG_INFINITY;
	double* M = new double[n];
	for (int v=0; v<n; v++) {
		M[v] = NEG_INFINITY;
		for (int k=0; k<n; k++) M[v] = min(M[v], (D[n][v]-D[k][v])/(n-k));

		lambda = max(lambda, M[v]);
	}

	return lambda;
}

/// /brief generate strongly connected componets for the digraph
/// The algorithm is implemented after "Cormen et al: Introduction to
/// agorithms", Chapter 22.5. It has a running time of O(V + E)
vector<Digraph*> GraphAlgorithms::generate_strongly_connected_components(Digraph* digraph) {
	vector<Digraph*> sccs;
	// step 1: call DFS(G) to compute finishing  times f[u] for each vertex u
	depth_first_search(digraph,digraph->node_vec);

	// step 2: compute G^T
	// step 3: call DFS(G^T). But in the main loop of DFS, 
	//         consider the verices in order of decreasing f[u] (as computed in step 1)
	// step 4: return the strongly connected components (trees in the depth-first formed in step 3)

	// bubble sort by decreasing order of f[u]
	int n = digraph->node_vec.size();
	Node** node_array = new Node*[n];
	int index = 0;
	for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) node_array[index++] = *iter;

	for (int i=1; i<n; i++) {
		for (int j=0; j<n-i; j++) {
			if(node_array[j]->f<node_array[j+1]->f) {
				Node* temp = node_array[j];
				node_array[j] = node_array[j+1];
				node_array[j+1] = temp;
			}
		}
	}
	// DFS on the transpose graph and generate the whole strongly connected components
	reverse_depth_first_search(digraph, node_array, n, sccs);


	// set the edges for these strongly connected components
	for (vector<Digraph*>::iterator sccs_iter = sccs.begin(); sccs_iter != sccs.end(); sccs_iter++) {
		Digraph* scc = *sccs_iter;
		for (vector<Node*>::iterator n_iter = scc->node_vec.begin(); n_iter != scc->node_vec.end(); n_iter++) {
			Node* node = *n_iter;

			for (list<Edge*>::iterator e_iter=node->out.begin(); e_iter != node->out.end(); e_iter++) {
				Edge* edge = *e_iter;
				Node* snk = edge->snk_node;
				if (find(scc->node_vec.begin(),scc->node_vec.end(),snk) != scc->node_vec.end()) 
					scc->add_scc_edge(edge);
			}
		}
	}

	delete[] node_array;
	return sccs;
}

void GraphAlgorithms::depth_first_search(Digraph* digraph, vector<Node*> new_node_vec) {
	if (new_node_vec.empty()) cout << "Empty vector. It should never arrive at here!" <<endl;

	vector<Node*>::iterator iter;
	for (iter = new_node_vec.begin(); iter != new_node_vec.end(); iter++) {
		Node* node = *iter;
		node->set_color(Node::WHITE);
		node->pi = NULL;
	}
	int time = 0;

	for (iter = new_node_vec.begin(); iter != new_node_vec.end(); iter++) {
		Node* node = *iter;
		if (node->get_color() == Node::WHITE)
			depth_first_search_visit(node,time);
	}
}

void GraphAlgorithms::depth_first_search_visit(Node* node, int &time){
	node->set_color(Node::GRAY);
	time++;
	node->d = time;
	for (list<Edge*>::iterator iter = node->out.begin(); iter != node->out.end(); iter++) {
		Edge* edge = *iter;
		Node* snk = edge->snk_node;
		if (snk->get_color() == Node::WHITE) {
			snk->pi = node;
			depth_first_search_visit(snk, time);
		}
	}
	node->set_color(Node::BLACK);
	time++;
	node->f = time;
}

void GraphAlgorithms::reverse_depth_first_search(Digraph* digraph, Node** node_array, int n, vector<Digraph*>& sccs) {
	if (n==0) cout << "Empty vector. It should never arrive at here!" <<endl;
	// We only need to identify color withou the discovery time and finish time
	for (int i=0; i<n; i++) {
		Node* node = node_array[i];
		node->set_color(Node::WHITE);
		node->pi = NULL;
	}

	for (int i=0; i<n; i++) {
		Node* node = node_array[i];
		if (node->get_color() == Node::WHITE) {
			Digraph* scc = new Digraph();
			reverse_depth_first_search_visit(node,scc);
			sccs.push_back(scc);
		}
	}
}

void GraphAlgorithms::reverse_depth_first_search_visit(Node* node, Digraph* scc){
	node->set_color(Node::GRAY);
	scc->add_node(node);
	for (list<Edge*>::iterator iter = node->in.begin(); iter != node->in.end(); iter++) {
		Edge* edge = *iter;
		Node* src = edge->src_node;
		if (src->get_color() == Node::WHITE) {
			src->pi = node;
			reverse_depth_first_search_visit(src,scc);
		}
	}
	node->set_color(Node::BLACK);
}

/// /brief generate strongly connected componets for the general digraph
/// The algorithm is implemented after "Cormen et al: Introduction to
/// agorithms", Chapter 22.5. It has a running time of O(V + E)
vector<GeneralDirectedGraph*> GraphAlgorithms::generate_strongly_connected_components(GeneralDirectedGraph* gDigraph) {
	vector<GeneralDirectedGraph*> sccs;
	// step 1: call DFS(G) to compute finishing  times f[u] for each vertex u
	depth_first_search(gDigraph,gDigraph->node_vec);

	// step 2: compute G^T
	// step 3: call DFS(G^T). But in the main loop of DFS, 
	//         consider the verices in order of decreasing f[u] (as computed in step 1)
	// step 4: return the strongly connected components (trees in the depth-first formed in step 3)

	// bubble sort by decreasing order of f[u]
	int n = gDigraph->node_vec.size();
	Node** node_array = new Node*[n];
	int index = 0;
	for (vector<Node*>::iterator iter = gDigraph->node_vec.begin(); iter != gDigraph->node_vec.end(); iter++) node_array[index++] = *iter;

	for (int i=1; i<n; i++) {
		for (int j=0; j<n-i; j++) {
			if(node_array[j]->f<node_array[j+1]->f) {
				Node* temp = node_array[j];
				node_array[j] = node_array[j+1];
				node_array[j+1] = temp;
			}
		}
	}
	// DFS on the transpose graph and generate the whole strongly connected components
	reverse_depth_first_search(gDigraph, node_array, n, sccs);


	// set the edges for these strongly connected components
	for (vector<GeneralDirectedGraph*>::iterator sccs_iter = sccs.begin(); sccs_iter != sccs.end(); sccs_iter++) {
		GeneralDirectedGraph* scc = *sccs_iter;
		for (vector<Node*>::iterator n_iter = scc->node_vec.begin(); n_iter != scc->node_vec.end(); n_iter++) {
			Node* node = *n_iter;

			for (list<Edge*>::iterator e_iter=node->out.begin(); e_iter != node->out.end(); e_iter++) {
				Edge* edge = *e_iter;
				Node* snk = edge->snk_node;
				if (find(scc->node_vec.begin(),scc->node_vec.end(),snk) != scc->node_vec.end()) 
					scc->add_scc_edge(edge);
			}
		}
	}
	delete[] node_array;
	return sccs;
}

void GraphAlgorithms::depth_first_search(GeneralDirectedGraph* gDigraph, vector<Node*> new_node_vec) {
	if (new_node_vec.empty()) cout << "Empty vector. It should never arrive at here!" <<endl;

	vector<Node*>::iterator iter;
	for (iter = new_node_vec.begin(); iter != new_node_vec.end(); iter++) {
		Node* node = *iter;
		node->set_color(Node::WHITE);
		node->pi = NULL;
	}
	int time = 0;

	for (iter = new_node_vec.begin(); iter != new_node_vec.end(); iter++) {
		Node* node = *iter;
		if (node->get_color() == Node::WHITE)
			depth_first_search_visit(node,time);
	}
}

void GraphAlgorithms::reverse_depth_first_search(GeneralDirectedGraph* digraph, Node** node_array, int n, vector<GeneralDirectedGraph*>& sccs) {
	if (n==0) cout << "Empty vector. It should never arrive at here!" <<endl;
	// We only need to identify color withou the discovery time and finish time
	for (int i=0; i<n; i++) {
		Node* node = node_array[i];
		node->set_color(Node::WHITE);
		node->pi = NULL;
	}

	for (int i=0; i<n; i++) {
		Node* node = node_array[i];
		if (node->get_color() == Node::WHITE) {
			GeneralDirectedGraph* scc = new GeneralDirectedGraph();
			reverse_depth_first_search_visit(node,scc);
			sccs.push_back(scc);
		}
	}
}

void GraphAlgorithms::reverse_depth_first_search_visit(Node* node, GeneralDirectedGraph* scc){
	node->set_color(Node::GRAY);
	scc->add_node(node);
	for (list<Edge*>::iterator iter = node->in.begin(); iter != node->in.end(); iter++) {
		Edge* edge = *iter;
		Node* src = edge->src_node;
		if (src->get_color() == Node::WHITE) {
			src->pi = node;
			reverse_depth_first_search_visit(src,scc);
		}
	}
	node->set_color(Node::BLACK);
}

double** GraphAlgorithms::generate_exec_request_matrix(Digraph* digraph) {
	int n = digraph->node_vec.size();
	double** A = new double*[n];
	for (int i=0; i<n; i++) A[i] = new double[n];

	for (int i=0; i<n; i++) for (int j=0; j<n; j++) A[i][j] = NEG_INFINITY;

	for (vector<Edge*>::iterator iter=digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
		Node* src = (*iter)->src_node;
		Node* snk = (*iter)->snk_node;

		int srcIndx = digraph->node_to_index[src];
		int snkIndx = digraph->node_to_index[snk];

		int weight = src->wcet;

		A[srcIndx][snkIndx] = weight;
	}

	for (vector<Node*>::iterator iter=digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
		int indx = digraph->node_to_index[*iter];
		double weight = 0.0;

		A[indx][indx] = max(A[indx][indx],weight);
	}
	
	return A;
}

double** GraphAlgorithms::generate_exec_request_matrix(Digraph* digraph,set<int> iSet) {
	int n = digraph->node_vec.size();
	double** A = new double*[n];
	for (int i=0; i<n; i++) A[i] = new double[n];

	for (int i=0; i<n; i++) for (int j=0; j<n; j++) A[i][j] = NEG_INFINITY;

	for (vector<Edge*>::iterator iter=digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
		Node* src = (*iter)->src_node;
		Node* snk = (*iter)->snk_node;

		int srcIndx = digraph->node_to_index[src];
		int snkIndx = digraph->node_to_index[snk];

		int weight = src->wcet;

		A[srcIndx][snkIndx] = weight;
	}

	for (set<int>::iterator iter=iSet.begin(); iter != iSet.end(); iter++) {
		int indx = *iter;
		double weight = 0.0;

		A[indx][indx] = max(A[indx][indx],weight);
	}
	
	return A;
}

void GraphAlgorithms::calculate_csum(Digraph* digraph) {
	int csum = 0;
	for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
		Node* node = *iter;
		csum += node->wcet;
	}
	digraph->c_sum = csum;
}

void GraphAlgorithms::calculate_tight_linear_bounds(Digraph* digraph) {
	int n = digraph->node_vec.size()*2;
	double** B = new double*[n];
	for (int i=0; i<n; i++) B[i] = new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) B[i][j] = NEG_INFINITY;

	vector<Node*>::iterator n_iter;
	list<Edge*>::iterator e_iter;
	for (n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		int index = digraph->node_to_index[*n_iter];
		B[index*2][index*2] = - digraph->linear_factor*digraph->pGCD;
	}

	for (n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		int srcIndx = digraph->node_to_index[node];
		for (e_iter = node->out.begin(); e_iter != node->out.end(); e_iter++) {
			Edge* edge = *e_iter;
			int snkIndx = digraph->node_to_index[edge->snk_node];
			//cout<<B[srcIndx*2][srcIndx*2+1]<<"\t"<<(double)node->wcet-digraph->linear_factor*digraph->gcd<<endl;
			B[srcIndx*2][srcIndx*2+1] = max(B[srcIndx*2][srcIndx*2+1] , (double)node->wcet-digraph->linear_factor*digraph->pGCD);
			//cout<<B[srcIndx*2][srcIndx*2+1]<<"\t"<<(double)node->wcet-digraph->linear_factor<<endl;
			B[srcIndx*2+1][snkIndx*2] = max(B[srcIndx*2+1][snkIndx*2], (1.0-(double)edge->separationTime/digraph->pGCD)*digraph->linear_factor*digraph->pGCD);
		}
	}

	//Utility::output_matrix(B,n,n);

	double** barB = MaxPlusAlgebra::calculate_metric_matrix(B,n);
	//Utility::output_matrix(barB,n,n);
	double W = 0;
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) W = max(W, barB[i][j]);
	
	double R = 0;
	for (int i=0; i<n; i+=2) for (int j=1; j<n; j+=2) {
		if (j%2==1) {
			double wcet;
			map<Node*, int>::iterator m_iter;
			for (m_iter = digraph->node_to_index.begin(); m_iter != digraph->node_to_index.end(); m_iter++) {
				if (m_iter->second == j/2) {
					wcet = m_iter->first->wcet;
					break;
				}
			}
			R = max(R, barB[i][j]-wcet*digraph->linear_factor);
		} 
		else
			R = max (R, barB[i][j]);
	}

	double X = 0;

	double a = NEG_INFINITY;
	double b = NEG_INFINITY;
	for (n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		for (e_iter = node->in.begin(); e_iter != node->in.end(); e_iter++) {
			Edge* edge = *e_iter;
			a = max(a, node->wcet-digraph->linear_factor*(edge->separationTime+node->deadline));
		}

		b = max(b, -digraph->linear_factor*node->deadline);
	}

	X = W+digraph->linear_factor*digraph->pGCD + min(a,b);

	digraph->c_rbf = W + digraph->linear_factor*digraph->pGCD;
	digraph->c_ibf = R + digraph->linear_factor*digraph->pGCD;
	digraph->c_dbf = X;

	// release
	for (int i=0; i<n; i++) delete[] B[i];
	delete[] B;

	for (int i=0; i<n; i++) delete[] barB[i];
	delete[] barB;
}

Digraph* GraphAlgorithms::generate_simple_digraph(Stateflow* sf) {
	Digraph* digraph = new Digraph(sf->scale);
	map<Transition*, Node*> tm;
	map<Node*,Transition*> n2t;

	// Generate nodes
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* tran = *iter;
		int index = sf->tran_index[tran];
		Node* node = new Node("simple_v"+Utility::int_to_string(index),sf->scale,tran->wcet,sf->gcd);
		digraph->add_node(node);
		tm[tran] = node;
		n2t[node] = tran;
	}

	// Generate edges
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* tran = *iter;
		Node* src = tm[tran];
		for (list<Transition*>::iterator iter1 = tran->snk->out.begin(); iter1 != tran->snk->out.end(); iter1++) {
			Transition* t = *iter1;
			Node* snk = tm[t];
			Edge* edge = new Edge(src, snk);
			edge->separationTime = Utility::math_gcd(tran->period, t->period);
			src->out.push_back(edge);
		    snk->in.push_back(edge);
			digraph->add_edge(edge);
		}
	}

	// add some null vertices for the vertices without out-edge
	bool containEmpty = false;
	Node* null_node = new Node("null_simple_v",sf->scale,0,0);
	for (vector<Node*>::iterator n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		if (!node->out.empty()) continue;
		containEmpty = true;		
		Edge* edge = new Edge(node, null_node);
		edge->separationTime = n2t[node]->period;
		node->out.push_back(edge);
		null_node->in.push_back(edge);
		digraph->add_edge(edge);
	}

	if (containEmpty) digraph->add_node(null_node);

	
	// set up the deadlines for all the vertices
	for (vector<Node*>::iterator n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		//int deadline = sf->hyperperiod; 
		if (node->out.empty()) continue;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator e_iter = node->out.begin(); e_iter != node->out.end(); e_iter++) {
			Edge* edge = *e_iter;
			deadline = min(deadline,edge->separationTime);
		}
		node->set_deadline(deadline);
	}
	
	tm.clear();
	n2t.clear();
	return digraph;
}

Digraph* GraphAlgorithms::generate_precise_digraph(Stateflow* sf) {
	Digraph* digraph = new Digraph(sf->scale);
	int n = sf->trans.size();
	int m = sf->n_time_instance-1;

	Node*** nodes = new Node**[n];
	for (int i=0; i<n; i++) nodes[i] = new Node*[m];

	// Generate nodes
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* tran = *iter;
		int i = sf->tran_index[tran];
		for (int j=0; j<m; j++) {
			int time = sf->index_time[j];
			if (time%tran->period == 0) {
				nodes[i][j] = new Node("precise_v"+Utility::int_to_string(i)+"_"+Utility::int_to_string(j), sf->scale,tran->wcet,sf->gcd);
			} else
				nodes[i][j] = NULL;
		}
	}

	// Generate edges
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			Node* src = nodes[j][i];
			if (src == NULL) continue;
			Transition* tran = sf->index_tran[j];
			for (list<Transition*>::iterator iter = tran->snk->out.begin(); iter != tran->snk->out.end(); iter++) {
				Transition* t = *iter;
				int index = sf->tran_index[t];
				bool isAfter = false; // used to identify whether the snk node is after the src node.
				for (int k=i+1; k<m; k++) {
					Node* snk = nodes[index][k];
					if (snk==NULL) continue;
					Edge* edge = new Edge(src,snk);
					edge->separationTime = sf->index_time[k]-sf->index_time[i];
					src->out.push_back(edge);
					snk->in.push_back(edge);
					digraph->add_edge(edge);
					isAfter = true;
					break; // find the first node
				}

				if (!isAfter) {
					for (int k=0; k<=i; k++) {
						Node* snk = nodes[index][k];
						if (snk==NULL) continue;
						Edge* edge = new Edge(src,snk);
						edge->separationTime = sf->index_time[m]-sf->index_time[i]+sf->index_time[k]-sf->index_time[0];
						src->out.push_back(edge);
						snk->in.push_back(edge);
						digraph->add_edge(edge);
						break; // find the first node
					}
				}
			}
		}
	}

	// add nodes
	int index = 0;
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			if (nodes[i][j] == NULL) continue;
			/*
			if (nodes[i][j]->in.empty() && nodes[i][j]->out.empty()) { 
				delete nodes[i][j]; 
				nodes[i][j] = NULL;
				continue;
			}
			*/
			nodes[i][j]->index = index++;
			digraph->add_node(nodes[i][j]);
		}
	}

	// add the null vertex for the vertices without out-edge
	bool containEmpty = false;
	Node* null_node = new Node("null_precise_v",sf->scale,0,0);
	for (vector<Node*>::iterator n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		if (!node->out.empty()) continue;
		containEmpty = true;
		Edge* edge = new Edge(node,null_node);
		int separationTime = 0;
		for (int i=0; i<n; i++) for (int j=0; j<m; j++) {
			if (nodes[i][j] == node) {
				separationTime = sf->trans.at(i)->period;
				break;
			}
		}
		edge->separationTime = separationTime;
		node->out.push_back(edge);
		null_node->in.push_back(edge);
		digraph->add_edge(edge);
	}
	if (containEmpty) digraph->add_node(null_node);

	// set up the deadlines for all the vertices
	for (vector<Node*>::iterator n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		if (node->out.empty()) continue;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator e_iter = node->out.begin(); e_iter != node->out.end(); e_iter++) {
			Edge* edge = *e_iter;
			deadline = min(deadline,edge->separationTime);
		}
		/*
		if (node->out.empty()) {
			for (int i=0; i<n; i++) for (int j=0; j<m; j++)
				if (nodes[i][j] == node) {
					deadline = sf->trans.at(i)->period;
					break;
				}
		}
		*/
		node->set_deadline(deadline);
	}

	for (int i=0; i<n; i++) delete[] nodes[i];
	delete[] nodes;

	return digraph;
}

bool GraphAlgorithms::isCyclic(Digraph* digraph) {
	// Mark all the vertices as not visited and not part of recursion
	// stack
	int n = digraph->node_vec.size();

	bool* visited = new bool[n];
	bool* recStack = new bool[n];

	for (int i=0; i<n; i++) {
		visited[i] = false;
		recStack[i] = false;
	}

	// Call the recursive helper function to detect cycle in different
	// DFS trees
	for (int i=0; i<n; i++)
		if (isCyclicUtil(digraph, i, visited, recStack))
			return true;

	return false;
}

bool GraphAlgorithms::isCyclicUtil(Digraph* digraph, int v, bool* visited, bool* recStack) {
	if (visited[v] == false) {
		// Mark the current node as visited and part of recursion stack
		visited[v] = true;
		recStack[v] = true;

		// Recur for all the vertices adjacent to this vertex
		Node* node = digraph->index_to_node[v]; 
		for (list<Edge*>::iterator iter = node->out.begin(); iter != node->out.end(); iter++) {
			Edge* edge = *iter;
			int w = digraph->node_to_index[edge->snk_node];
			if (!visited[w] && isCyclicUtil(digraph,w,visited,recStack))
				return true;
			else if (recStack[w])
				return true;
		}
	}

	recStack[v] = false; // remove the vertex from recursion stack
	return false;
}

bool GraphAlgorithms::isCyclic(Stateflow* sf) {
	// Mark all the vertices as not visited and not part of recursion
	// stack
	int n = sf->states.size();
	bool* visited = new bool[n];
	bool* recStack = new bool[n];

	for (int i=0; i<n; i++) {
		visited[i] = false;
		recStack[i] = false;
	}

	// Call the recursive helper function to detect cycle in different
	// DFS trees
	for (int i=0; i<n; i++)
		if (isCyclicUtil(sf, i, visited, recStack)) {
			// release visisted and recStack
			delete[] visited;
			delete[] recStack;

			return true;
		}

	delete[] visited;
	delete[] recStack;
	return false;
}

bool GraphAlgorithms::isCyclicUtil(Stateflow* sf, int v, bool* visited, bool* recStack) {
	if (visited[v] == false) {
		// Mark the current node as visited and part of recursion stack
		visited[v] = true;
		recStack[v] = true;

		// Recur for all the vertices adjacent to this vertex
		State* state = sf->index_state[v]; 
		for (list<Transition*>::iterator iter = state->out.begin(); iter != state->out.end(); iter++) {
			Transition* t = *iter;
			int w = sf->state_index[t->snk];
			if (!visited[w] && isCyclicUtil(sf,w,visited,recStack))
				return true;
			else if (recStack[w])
				return true;
		}
	}

	recStack[v] = false; // remove the vertex from recursion stack
	return false;
}