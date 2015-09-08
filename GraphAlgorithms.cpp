#include "GraphAlgorithms.h"

#include "MaxPlusAlgebra.h"
#include <algorithm>

/// /brief an implementation of the Karp's algorithm
/// R.M. Karp, A characterization of the minimum cycle mean in a digraph, Discrete Math. 23 (1978) 309�11.
/// Note: this function is just right when operating on UDRT in case of the length.
double GraphAlgorithms::calculate_maximum_cycle_mean(Digraph* digraph) {
	digraph->calculate_gcd();
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
	lambda = lambda/digraph->gcd;
	digraph->linear_factor = lambda;
	return lambda;
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
		B[index*2][index*2] = - digraph->linear_factor;
	}

	for (n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		int srcIndx = digraph->node_to_index[node];
		for (e_iter = node->out.begin(); e_iter != node->out.end(); e_iter++) {
			Edge* edge = *e_iter;
			int snkIndx = digraph->node_to_index[edge->snk_node];
			B[srcIndx*2][srcIndx*2+1] = max(B[srcIndx*2][srcIndx*2+1] , (double)node->wcet-digraph->linear_factor);
			B[srcIndx*2+1][snkIndx*2] = max(B[srcIndx*2+1][snkIndx*2], (1.0-(double)edge->separationTime/digraph->gcd)*digraph->linear_factor);
		}
	}

	double** barB = MaxPlusAlgebra::calculate_metric_matrix(B,n);
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
			R = max(R, barB[i][j]-wcet*digraph->linear_factor/digraph->gcd);
		} 
		else
			R = max (R, barB[i][j]);
	}

	double X = 0;

	double a = NEG_INFINITY;
	double b = NEG_INFINITY;
	for (n_iter = digraph->node_vec.begin(); n_iter != digraph->node_vec.end(); n_iter++) {
		Node* node = *n_iter;
		int srcIndx = digraph->node_to_index[node];
		for (e_iter = node->out.begin(); e_iter != node->out.end(); e_iter++) {
			Edge* edge = *e_iter;
			a = max(a, node->wcet-digraph->linear_factor*(edge->separationTime+node->deadline)/digraph->gcd);
		}

		b = max(b, -digraph->linear_factor*node->deadline/digraph->gcd);
	}

	X = W+digraph->linear_factor + max(a,b);

	digraph->c_rbf = W + digraph->linear_factor;
	digraph->c_ibf = R + digraph->linear_factor;
	digraph->c_dbf = X;
}

Digraph* GraphAlgorithms::generate_simple_digraph(Stateflow* sf) {
	Digraph* digraph = new Digraph(sf->scale);
	map<Transition*, Node*> tm;

	// Generate nodes
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* tran = *iter;
		int index = sf->tran_index[tran];
		Node* node = new Node("simple_v"+Utility::int_to_string(index),sf->scale,tran->wcet,tran->period);
		digraph->add_node(node);
		tm[tran] = node;
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
				nodes[i][j] = new Node("precise_v"+Utility::int_to_string(i)+"_"+Utility::int_to_string(j), sf->scale,tran->wcet,tran->period);
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
	for (int i=0; i<n; i++)
		for (int j=0; j<m; j++) {
			if (nodes[i][j] == NULL) continue;
			if (nodes[i][j]->in.empty() && nodes[i][j]->out.empty()) continue;
			nodes[i][j]->index = index++;
			digraph->add_node(nodes[i][j]);
		}

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
		if (isCyclicUtil(sf, i, visited, recStack))
			return true;

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