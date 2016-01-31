/* \file Digraph.cpp
*  this file define method functions of digraph, unitdigraph and grandigraph class. 
*  \author Chao Peng
*  
*  changes
*  ------
*  21-aug-2015 : initial revision (CP)
*
*/

#include "Digraph.h"

#include "MaxPlusAlgebra.h"
#include "GraphAlgorithms.h"

extern double EPSILON;

using namespace std;

// ===========================================================================
// Method functions of Digraph class
// ===========================================================================
Digraph::~Digraph() {
	//cout<<"Digraph Destruction Start"<<endl;
	// release sccs
	for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
		(*iter)->node_vec.clear();
		(*iter)->edge_vec.clear();
		delete *iter;
		*iter = NULL;
	}
	sccs.clear();
	
	// release nodes
	//cout<<node_vec.empty()<<endl;
	if (iNode != 0) {
		for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
	}
	node_vec.clear();
	cnode_vec.clear();

	// release edges
	if (iEdge!=0) {
		for (vector<Edge*>::iterator iter = edge_vec.begin(); iter != edge_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
	}
	edge_vec.clear();

	node_to_index.clear();
	index_to_node.clear();
	edge_to_index.clear();
	index_to_edge.clear();

	if(unit_digraph != NULL) delete unit_digraph;
	if(gran_digraph != NULL) delete gran_digraph;

	//cout<<"Digraph End"<<endl;
}

void Digraph::prepare_digraph() {
	 generate_strongly_connected_components();
	 check_strongly_connected();
	 calculate_period_gcd();
	 calculate_all_gcd();
	 calculate_linear_factor();
	 calculate_linear_upper_bounds();
}

void Digraph::generate_strongly_connected_components() {
	this->sccs = GraphAlgorithms::generate_strongly_connected_components(this);
}

/// /brief generate highly connected components from the SCCS of this digraph
void Digraph::generate_highly_connected_components() {
	vector<Digraph*>::iterator iter;
	for (iter = sccs.begin(); iter != sccs.end(); iter++) {
		Digraph* scc = *iter;
		if (scc->linear_factor == this->linear_factor)
			this->hccs.push_back(scc);
	}
}

void Digraph::check_strongly_connected() {
	if (sccs.size() == 1) this->strongly_connected = true;
	else this->strongly_connected = false;
}

void Digraph::calculate_linear_factor() {
	//this->generate_strongly_connected_components();
	double lambda = 0.0;
	vector<Digraph*>::iterator iter;
	for (iter = sccs.begin(); iter != sccs.end(); iter++) {
		GraphAlgorithms::calculate_maximum_cycle_mean(*iter);
		lambda = max(lambda,(*iter)->linear_factor);
	}
	this->linear_factor = lambda;
}

void Digraph::calculate_linear_factor2() {
	this->linear_factor = GraphAlgorithms::calculate_maximum_cycle_mean2(this);
}

void Digraph::calculate_csum() {
	// calculate C^sum
	GraphAlgorithms::calculate_csum(this);
}

void Digraph::calculate_linear_upper_bounds() {
	// calculate C^rbf, C^ibf and C^dbf
	GraphAlgorithms::calculate_tight_linear_bounds(this);
}

void Digraph::scale_wcet(double factor) {
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		Node* node = *iter;
		int wcet = (int)(factor*node->wcet);
		if (wcet == 0) wcet = rand()%5+1; // force the wcet not to be 0
		node->wcet = wcet;
	}
}


void Digraph::calculate_tf0(Digraph** digraphs, int i) {
	double util = 0;
	double sum = 0;
	for (int k=0; k<=i; k++) {
		Digraph* digraph = digraphs[k];
		util += digraph->linear_factor;

		sum += digraph->c_sum*2;
		if (k==i) sum -= digraph->c_sum;
	}

	if (util >= 1) tf0 = POS_INFINITY;

	tf0 = sum/(1-util);
}

void Digraph::calculate_tf1(Digraph** digraphs, int i) {
	double util = digraphs[i]->linear_factor;
	double sum = digraphs[i]->c_dbf;
	for (int k=0; k<i; k++) {
		Digraph* digraph = digraphs[k];
		util += digraph->linear_factor;

		sum += digraph->c_rbf;
	}

	if (util >= 1) tf1 = POS_INFINITY;

	tf1 = sum/(1-util);
}

void Digraph::calculate_tf2(Digraph** digraphs, int i) {
	double util = digraphs[i]->linear_factor;
	double sum = digraphs[i]->c_dbf;
	for (int k=0; k<i; k++) {
		Digraph* digraph = digraphs[k];
		util += digraph->linear_factor;

		sum += digraph->c_ibf;
	}

	if (util >= 1) tf2 = POS_INFINITY;

	tf2 = sum/(1-util);
}

void Digraph::prepare_rbf_calculation(bool debug) {
	unit_digraph = new UnitDigraph(this);
	unit_digraph->prepare(debug);
}

void Digraph::prepare_rbf_calculation_without_periodicity(bool debug) {
	unit_digraph = new UnitDigraph(this);
	unit_digraph->prepare_without_periodicity(debug);
}

void Digraph::prepare_ibf_calculation(bool debug) {
	gran_digraph = new GranDigraph(this);
	gran_digraph->prepare(debug);
}

double Digraph::rbf(int t) {
	if (t==0) return 0;
	double rbf = unit_digraph->get_rbf(t);
	double rbf_leaf = 0; // if some nodes have no out edge
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++)
		rbf_leaf = max(rbf_leaf, (double)(*iter)->wcet);

	return max(rbf,rbf_leaf);
}

double** Digraph::rbf_exec_req_matrix(int t, int& n) {
	n = unit_digraph->n_size;
	return unit_digraph->matrix_map[t];
}

double Digraph::ibf(int t) {
	return gran_digraph->get_ibf(t);
}

double** Digraph::ibf_exec_req_matrix(int t, int& n) {
	n = gran_digraph->n_size;
	return gran_digraph->matrix_map[t];
}

double Digraph::dbf(int t) {
	int large_t = 100; // sufficiently large t
	double factor = (double) t/pGCD;
	t = floor(factor);
	//write_graphviz(cout);
	if (dbf_map.find(t) != dbf_map.end())
		return dbf_map[t];
	if (t < large_t)
		return calculate_dbf(t);
	else 
		return unit_digraph->get_dbf(t);
}

double Digraph::calculate_dbf(int t) {
	vector<vector<DemandTriple*>> DT;
	// initial DT0
	vector<DemandTriple*> DT0;
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		Node* node = *iter;
		if (node->deadline == INT_MAX) continue;
		DemandTriple* dt = new DemandTriple(node->wcet, node->deadline, node);
		DT0.push_back(dt);
	}
	DT.push_back(DT0);

	for (int k=1; k<=t; k++) {
		vector<DemandTriple*> DTk;
		vector<DemandTriple*> PrevDT = DT.back();
		for (vector<DemandTriple*>::iterator iter = PrevDT.begin(); iter != PrevDT.end(); iter++) {
			DemandTriple* dt = *iter;
			Node* src = dt->node;
			for (list<Edge*>::iterator eIter = src->out.begin(); eIter != src->out.end(); eIter++) {
				Edge* edge = *eIter;
				Node* snk = edge->snk_node;
				if (snk->deadline == INT_MAX) continue;
				int e = dt->e + snk->wcet;
				int d = dt->d - src->deadline + edge->separationTime + snk->deadline;

				bool dominated = false;
				for (vector<DemandTriple*>::iterator it = DTk.begin(); it != DTk.end(); it++) {
					DemandTriple* dt2 = *it;
					if (dt2->e >= e && dt2->d <= d && snk == dt2->node) {
						dominated = true;
						break;
					}
					if (dt2->e <= e && dt2->d >= d && snk == dt2->node) {
						dt2->e = e;
						dt2->d = d;
						dominated = true;
						break;
					}
				}

				if (!dominated && d<=t*pGCD) {
					DemandTriple* ndt = new DemandTriple(e,d,snk);
					DTk.push_back(ndt);
				}
			}
		}
		if (!DTk.empty())
			DT.push_back(DTk);
	}

	//DT.erase(DT.begin()); // remove thd DT0
	int max = 0;
	for (vector<vector<DemandTriple*>>::iterator iter = DT.begin(); iter != DT.end(); iter++) {
		vector<DemandTriple*> DTk = *iter;
		for (vector<DemandTriple*>::iterator it = DTk.begin(); it != DTk.end(); it++) {
			DemandTriple* dt = *it;
			if (dt->d <= t*pGCD)
				max = std::max(max,dt->e);
		}
	}
	dbf_map[t] = max;

	// delete DemandTriple
	for (vector<vector<DemandTriple*>>::iterator iter = DT.begin(); iter != DT.end(); iter++) {
		vector<DemandTriple*> DTk = *iter;
		for (vector<DemandTriple*>::iterator it = DTk.begin(); it != DTk.end(); it++) {
			delete *it;
			*it = NULL;
		}
		DTk.clear();
	}
	DT.clear();

	return max;
}

void Digraph::write_graphviz(ostream& out) {
	out << "digraph G {" <<endl;
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		Node* node = *iter;
		out << node->name << " [label=\" " << node->name << " / " << node->wcet << " / " << node->deadline << " \"]" << endl;
	}
	for (vector<Edge*>::iterator iter = edge_vec.begin(); iter != edge_vec.end(); iter++) {
		Edge* edge = *iter;
		out << edge->src_node->name << " -> " << edge->snk_node->name 
			<< " [label=\" " << edge->separationTime << " \"]" << endl;
	}
	out << "}" <<endl;
}

// ===========================================================================
// Method functions of UnitDigraph class
// ===========================================================================
UnitDigraph::~UnitDigraph() {
	//cout<<"UnitDigraph Destruction Start"<<endl;

	// release nodes
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	node_vec.clear();

	// release edges
	for (vector<Edge*>::iterator iter = edge_vec.begin(); iter != edge_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	edge_vec.clear();

	node_to_index.clear();
	index_to_node.clear();
	edge_to_index.clear();
	index_to_edge.clear();

	// release sccs
	for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	sccs.clear();

	if (unit_digraph != NULL) delete unit_digraph;
	if (gran_digraph != NULL) delete gran_digraph;

	iSet.clear();
	original_node_to_index.clear();
	index_to_original_node.clear();

	// release matrix
	/*
	for (int i=0; i<n_size; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
	*/

	// delete execution request matrix
	for (map<int, double**>::iterator iter = matrix_map.begin(); iter != matrix_map.end(); iter++) {
		double** temp = iter->second;
		for (int i=0; i<n_size; i++)
			delete[] temp[i];
		delete[] temp;
	}

	matrix_map.clear();
	maximum_element_map.clear();

	//cout<<"UnitDigraph End"<<endl;
}

/// \brief prepare function
void UnitDigraph::prepare(bool debug) {
	// calculate gcd
	//calculate_period_gcd();
	//gcd = pGCD;
	// generate UDRT
	generate_unit_digraph();
	n_size = node_vec.size();

	// generate execution request matrix
	generate_exec_request_matrix();

	// Calculate linear period
	calculate_linear_period(debug);

	// Calculate linear defect
	calculate_linear_defect();

	// Calculate the execution request matrix power
	calculate_exec_request_matrix_power(tf);
}

void UnitDigraph::prepare_without_periodicity(bool debug) {
	generate_unit_digraph();
	n_size = node_vec.size();

	// generate execution request matrix
	generate_exec_request_matrix();

	// Calculate the execution request matrix power
	calculate_exec_request_matrix_power_without_periodicity(tf);
}

void UnitDigraph::generate_unit_digraph() {
	for (vector<Node*>::iterator iter = origin->node_vec.begin(); iter != origin->node_vec.end(); iter++) {
		Node* node = *iter;

		int maxSeparationtime = this->gcd;
		for (list<Edge*>::iterator iter = node->out.begin(); iter != node->out.end(); iter++) {
			Edge* edge = *iter;
			maxSeparationtime = std::max(maxSeparationtime, edge->separationTime);
		}

		int numNode = maxSeparationtime/this->gcd;
		node->unitNodes = new Node*[numNode];
		node->unitNodeNum = numNode;

		// Create new nodes
		for (int j=0; j<numNode; j++) {
			string name = "u_"+node->name+"_"+Utility::int_to_string(j);
			int deadline = 1; // not important
			int wcet = 0;
			if (j==0) {
				wcet = node->wcet;
				deadline = node->deadline; // the demand bound functions of the original digraph and UDRT are same
			}

			Node* newNode = new Node(name,node->scale,wcet,deadline);

			this->add_node(newNode, node);

			node->unitNodes[j] = newNode;
		}

		// Create new edges
		for (int j=1; j<numNode; j++) {
			Node* src = node->unitNodes[j-1];
			Node* snk = node->unitNodes[j];
			Edge* newEdge = new Edge(src, snk);
			newEdge->set_separation_time(1);
			src->out.push_back(newEdge);
			snk->in.push_back(newEdge);

			this->add_edge(newEdge);
		}
	}

	for (int i=0; i<origin->edge_vec.size(); i++) {
		Edge* edge = origin->edge_vec.at(i);
		int separation = edge->separationTime/this->gcd;
		Node* src = edge->src_node->unitNodes[separation-1];
		Node* snk = edge->snk_node->unitNodes[0];

		Edge* newEdge = new Edge(src, snk);
		newEdge->set_separation_time(1);
		src->out.push_back(newEdge);
		snk->in.push_back(newEdge);

		this->add_edge(newEdge);
	}
}

void UnitDigraph::scc_generate_unit_digraph() {
	for (vector<Node*>::iterator iter = origin->node_vec.begin(); iter != origin->node_vec.end(); iter++) {
		Node* node = *iter;

		int maxSeparationtime = this->gcd;
		for (list<Edge*>::iterator iter = node->scc_out.begin(); iter != node->scc_out.end(); iter++) {
			Edge* edge = *iter;
			maxSeparationtime = std::max(maxSeparationtime, edge->separationTime);
		}

		int numNode = maxSeparationtime/this->gcd;
		node->unitNodes = new Node*[numNode];
		node->unitNodeNum = numNode;

		// Create new nodes
		for (int j=0; j<numNode; j++) {
			string name = "u_"+node->name+"_"+Utility::int_to_string(j);
			int deadline = 1; // not important
			int wcet = 0;
			if (j==0) {
				wcet = node->wcet;
				deadline = node->deadline; // the demand bound functions of the original digraph and UDRT are same
			}

			Node* newNode = new Node(name,node->scale,wcet,deadline);

			this->add_node(newNode, node);

			node->unitNodes[j] = newNode;
		}

		// Create new edges
		for (int j=1; j<numNode; j++) {
			Node* src = node->unitNodes[j-1];
			Node* snk = node->unitNodes[j];
			Edge* newEdge = new Edge(src, snk);
			newEdge->set_separation_time(1);
			src->out.push_back(newEdge);
			snk->in.push_back(newEdge);

			this->add_edge(newEdge);
		}
	}

	for (int i=0; i<origin->edge_vec.size(); i++) {
		Edge* edge = origin->edge_vec.at(i);
		int separation = edge->separationTime/this->gcd;
		Node* src = edge->src_node->unitNodes[separation-1];
		Node* snk = edge->snk_node->unitNodes[0];

		Edge* newEdge = new Edge(src, snk);
		newEdge->set_separation_time(1);
		src->out.push_back(newEdge);
		snk->in.push_back(newEdge);

		this->add_edge(newEdge);
	}

	for (vector<Node*>::iterator iter = origin->node_vec.begin(); iter != origin->node_vec.end(); iter++) {
		Node* node = *iter;

		delete[] node->unitNodes;
		node->unitNodes = NULL;
	}
}

void UnitDigraph::generate_exec_request_matrix() {
	matrix = GraphAlgorithms::generate_exec_request_matrix(this);
	matrix_map[1]= matrix; // the execution request matrix of (0,gcd]
	double maximum_element = MaxPlusAlgebra::maximum_element(matrix, n_size);
	int temp = static_cast<int> (maximum_element);
	maximum_element_map[1] = temp;
	//matrix_map.insert(pair<int, double**>(1,matrix)); 
}

void UnitDigraph::calculate_linear_period(bool debug) {
	if (origin->strongly_connected)
		lper = MaxPlusAlgebra::calculate_linear_period(matrix,n_size,lfac*gcd,debug);
}

void UnitDigraph::calculate_linear_defect() {
	if (origin->strongly_connected)
		ldef = MaxPlusAlgebra::calculate_linear_defect(matrix_map, maximum_element_map, n_size, lfac, lper, tf, gcd);
}

void UnitDigraph::calculate_exec_request_matrix_power(int tf) {
	/// Calculate all the execution request matrix from 2 to tf
	/*
	for (int t=2; t<=tf; t++) {
		if (matrix_map.find(t) != matrix_map.end())
			continue;
		double** matrix_power;
		if (t>ldef+lper && origin->strongly_connected)
			matrix_power = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(matrix_map[t-lper],n_size,lfac*gcd,lper);
		else
			matrix_power = MaxPlusAlgebra::multiply_maxplus_matrix(matrix_map[t-1],matrix_map[1],n_size);
		matrix_map[t] = matrix_power;
		double maximum_element = MaxPlusAlgebra::maximum_element(matrix_power, n_size);
		int temp = static_cast<int> (maximum_element);
		maximum_element_map[t] = temp;
	}
	*/

	/// Calculate all the execution request matrix from 2 to tf
	/// Use rbf(t+p)= rbf(t)+p*q to calculate rbf(t) instead of the matrix operaition
	for (int t=2; t<=tf; t++) {
		if (matrix_map.find(t) != matrix_map.end())
			continue;

		if (t>ldef+lper && origin->strongly_connected && maximum_element_map.find(t-lper)!=maximum_element_map.end()) {
			double maximum_element = maximum_element_map[t-lper]+lfac*gcd*lper;
			int temp = static_cast<int> (maximum_element);
			maximum_element_map[t] = temp;
			continue;
		}

		double** matrix_power;
		matrix_power = MaxPlusAlgebra::multiply_maxplus_matrix(matrix_map[t-1],matrix_map[1],n_size);
		matrix_map[t] = matrix_power;
		double maximum_element = MaxPlusAlgebra::maximum_element(matrix_power, n_size);
		int temp = static_cast<int> (maximum_element);
		maximum_element_map[t] = temp;
	}

}

void UnitDigraph::calculate_exec_request_matrix_power_without_periodicity(int tf) {
	for (int t=2; t<=tf; t++) {
		if (matrix_map.find(t) != matrix_map.end())
			continue;
		double** matrix_power = MaxPlusAlgebra::multiply_maxplus_matrix(matrix_map[t-1],matrix_map[1],n_size);
		matrix_map[t] = matrix_power;
		double maximum_element = MaxPlusAlgebra::maximum_element(matrix_power, n_size);
		int temp = static_cast<int> (maximum_element);
		maximum_element_map[t] = temp;
	}
}

int UnitDigraph::get_rbf(int t) {
	// we should consider the gcd
	double factor = (double) t/gcd;
	t = ceil(factor);
	if (maximum_element_map.find(t) != maximum_element_map.end())
		return maximum_element_map[t];
	else 
		//calculate_exec_request_matrix_power_without_periodicity(t);
		calculate_exec_request_matrix_power(t);
	return maximum_element_map[t];
}

int UnitDigraph::calculate_demand_bound_function(int i, int j, int t) {
	Node* first = index_to_original_node[i];
	Node* last = index_to_original_node[j];
	int tprim = 1000*scale;
	for (list<Edge*>::iterator iter = last->in.begin(); iter != last->in.end(); iter++) {
		Edge* edge = *iter;
		tprim = min(tprim,edge->separationTime/gcd);
	}
	t = t - last->deadline/gcd - tprim;
	if (t<0) return 0;
	if (t==0) return last->wcet;

	double** t_matrix = matrix_map[t];
	return t_matrix[i][j]+last->wcet;
}

int UnitDigraph::get_dbf(int t) {
	int max = 0;
	for (set<int>::iterator i_iter=iSet.begin(); i_iter!= iSet.end(); i_iter++) {
		int i = *i_iter;
		for (set<int>::iterator j_iter=iSet.begin(); j_iter!= iSet.end(); j_iter++) {
			int j = *j_iter;
			max = std::max(max, calculate_demand_bound_function(i,j,t));
		}
	}
	return max;
}

// ===========================================================================
// Method functions of GranDigraph class
// ===========================================================================
GranDigraph::~GranDigraph() {
	//cout<<"GranDigraph Destruction Start"<<endl;

	// release nodes
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	node_vec.clear();

	// release edges
	for (vector<Edge*>::iterator iter = edge_vec.begin(); iter != edge_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	edge_vec.clear();

	node_to_index.clear();
	index_to_node.clear();
	edge_to_index.clear();
	index_to_edge.clear();

	// release sccs
	for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	sccs.clear();

	if (unit_digraph != NULL) delete unit_digraph;
	if (gran_digraph != NULL) delete gran_digraph;

	// release matrix
	/*
	for (int i=0; i<n_size; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
	*/

	// release matrix power
	// TODO:

	matrix_map.clear();
	maximum_element_map.clear();

	//cout<<"GranDigraph End"<<endl;
}

/// \brief prepare function
void GranDigraph::prepare(bool debug) {

	// calculate gcd
	//calculate_all_gcd();
	//gcd = aGCD;
	// generate GDRT
	generate_gran_digraph();
	n_size = node_vec.size();

	// generate execution request matrix
	generate_exec_request_matrix();

	// Calculate linear period
	calculate_linear_period(debug);

	// Calculate linear defect
	calculate_linear_defect();

	// Calculate the execution request matrix power
	calculate_exec_request_matrix_power(tf);
}

void GranDigraph::generate_gran_digraph() {
	for (vector<Node*>::iterator iter = origin->node_vec.begin(); iter != origin->node_vec.end(); iter++) {
		Node* node = *iter;

		int maxSeparationtime = this->gcd;
		for (list<Edge*>::iterator iter = node->out.begin(); iter != node->out.end(); iter++) {
			Edge* edge = *iter;
			maxSeparationtime = std::max(maxSeparationtime, edge->separationTime);
		}

		int numNode = maxSeparationtime/this->gcd;
		node->unitNodes = new Node*[numNode];

		// Create new nodes
		for (int j=0; j<numNode; j++) {
			string name = "g_"+node->name+"_"+Utility::int_to_string(j);
			int deadline = 1; // not important
			int wcet = 0;
			if (j<node->wcet/this->gcd) wcet = this->gcd;

			Node* newNode = new Node(name,node->scale,wcet,deadline);

			this->add_node(newNode);

			node->unitNodes[j] = newNode;
		}

		// Create new edges
		for (int j=1; j<numNode; j++) {
			Node* src = node->unitNodes[j-1];
			Node* snk = node->unitNodes[j];
			Edge* newEdge = new Edge(src, snk);
			newEdge->set_separation_time(1);
			src->out.push_back(newEdge);
			snk->in.push_back(newEdge);

			this->add_edge(newEdge);
		}
	}

	for (int i=0; i<origin->edge_vec.size(); i++) {
		Edge* edge = origin->edge_vec.at(i);
		int separation = edge->separationTime/this->gcd;
		Node* src = edge->src_node->unitNodes[separation-1];
		Node* snk = edge->snk_node->unitNodes[0];

		Edge* newEdge = new Edge(src, snk);
		newEdge->set_separation_time(1);
		src->out.push_back(newEdge);
		snk->in.push_back(newEdge);

		this->add_edge(newEdge);
	}
}

void GranDigraph::generate_exec_request_matrix() {
	matrix = GraphAlgorithms::generate_exec_request_matrix(this);
	matrix_map[1]= matrix; // the execution request matrix of (0,gcd]
	double maximum_element = MaxPlusAlgebra::maximum_element(matrix, n_size);
	int temp = static_cast<int> (maximum_element);
	maximum_element_map[1] = temp;
	//matrix_map.insert(pair<int, double**>(1,matrix)); 
}

void GranDigraph::calculate_linear_period(bool debug) {
	lper = MaxPlusAlgebra::calculate_linear_period(matrix,n_size,lfac, debug);
}

void GranDigraph::calculate_linear_defect() {
	ldef = MaxPlusAlgebra::calculate_linear_defect(matrix_map, maximum_element_map, n_size, lfac, lper, tf,gcd);
}

void GranDigraph::calculate_exec_request_matrix_power(int tf) {
	/// Calculate all the execution request matrix from 2 to tf
	for (int t=2; t<=tf; t++) {
		if (matrix_map.find(t) != matrix_map.end())
			continue;
		double** matrix_power;
		if (t>ldef+lper)
			matrix_power = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(matrix_map[t-lper],n_size,lfac*gcd,lper);
		else
			matrix_power = MaxPlusAlgebra::multiply_maxplus_matrix(matrix_map[t-1],matrix_map[1],n_size);
		matrix_map[t] = matrix_power;
		double maximum_element = MaxPlusAlgebra::maximum_element(matrix_power, n_size);
		int temp = static_cast<int> (maximum_element);
		maximum_element_map[t] = temp;
	}
}

int GranDigraph::get_ibf(int t) {
	if (t==0) return 0;
	// we should consider the gcd
	double factor = (double) t/gcd;
	int t_gcd = ceil(factor);
	int a = maximum_element_map[t_gcd];
	int b;
	if (t_gcd-1==0) b=0; else b = maximum_element_map[t_gcd-1];
	if (abs(a-b)<EPSILON) 
		return b;
	else {
		//std::cout << "test:" <<(a == b+gcd) <<std::endl;
		return b+t-(t_gcd-1)*gcd;
	}
}


	