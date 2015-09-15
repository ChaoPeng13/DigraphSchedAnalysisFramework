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
	 calculate_gcd();
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

void Digraph::calculate_linear_upper_bounds() {
	// calculate C^sum
	GraphAlgorithms::calculate_csum(this);
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

void Digraph::prepare_ibf_calculation(bool debug) {
	gran_digraph = new GranDigraph(this);
	gran_digraph->prepare(debug);
}

double Digraph::rbf(int t) {
	return unit_digraph->get_rbf(t);
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
	double factor = (double) t/gcd;
	t = floor(factor);
	return unit_digraph->get_dbf(t);
}

void Digraph::write_graphviz(ostream& out) {
	out << "digraph G {" <<endl;
	for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
		Node* node = *iter;
		out << node->name << " [label=\"" << node->name << "/" << node->wcet << "\"]" << endl;
	}
	for (vector<Edge*>::iterator iter = edge_vec.begin(); iter != edge_vec.end(); iter++) {
		Edge* edge = *iter;
		out << edge->src_node->name << " -> " << edge->snk_node->name 
			<< " [label=\"" << edge->separationTime << "\"]" << endl;
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

	// release matrix power
	// TODO:

	matrix_map.clear();
	maximum_element_map.clear();

	//cout<<"UnitDigraph End"<<endl;
}

/// \brief prepare function
void UnitDigraph::prepare(bool debug) {
	// calculate gcd
	calculate_gcd();

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
		ldef = MaxPlusAlgebra::calculate_linear_defect(matrix_map, maximum_element_map, n_size, lfac*gcd, lper*gcd, origin->tf);
}

void UnitDigraph::calculate_exec_request_matrix_power(int tf) {
	/// Calculate all the execution request matrix from 2 to tf
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
}

int UnitDigraph::get_rbf(int t) {
	// we should consider the gcd
	double factor = (double) t/gcd;
	t = ceil(factor);
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
	calculate_gcd();
	
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
	ldef = MaxPlusAlgebra::calculate_linear_defect(matrix_map, maximum_element_map, n_size, lfac, lper, origin->tf);
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


	