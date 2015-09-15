#if 0
#include <gtest/gtest.h>

#include "Digraph.h"
#include "MaxPlusAlgebra.h"

Digraph* digraph;

#define output false
extern double EPSILON;

/**
 * Used to test deep first search
 * The digraph comes from Figure 22-4 of "Cormen et al: Introduction to
 * agorithms", Chapter 22.3.
 */
void generateDigraph0() {
	digraph = new Digraph();

	// creat Node
	Node* u = new Node("u", 1, 1,1,4);
	Node* v = new Node("v", 2, 1,2,4);
	Node* w = new Node("w", 3, 1,2,4);
	Node* x = new Node("x", 4, 1,3,4);
	Node* y = new Node("y", 5, 1,2,4);
	Node* z = new Node("z", 6, 1,3,4);

	digraph->add_node(u);
	digraph->add_node(v);
	digraph->add_node(w);
	digraph->add_node(x);
	digraph->add_node(y);
	digraph->add_node(z);

	// Create Edge
	Edge* uv = new Edge(u,v);
	Edge* ux = new Edge(u,x);
	Edge* vy = new Edge(v,y);
	Edge* wy = new Edge(w,y);
	Edge* wz = new Edge(w,z);
	Edge* xv = new Edge(x,v);
	Edge* yx = new Edge(y,x);
	Edge* zz = new Edge(z,z);

	uv->set_separation_time(4);
	ux->set_separation_time(4);
	vy->set_separation_time(4);
	wy->set_separation_time(4);
	wz->set_separation_time(4);
	xv->set_separation_time(4);
	yx->set_separation_time(4);
	zz->set_separation_time(4);

	digraph->add_edge(uv);
	digraph->add_edge(ux);
	digraph->add_edge(vy);
	digraph->add_edge(wy);
	digraph->add_edge(wz);
	digraph->add_edge(xv);
	digraph->add_edge(yx);
	digraph->add_edge(zz);
	
	if (output) {
		digraph->write_graphviz(std::cout);
	}
}


/**
 * Used to test strongly connected components algortihm
 * The digraph comes from Figure 22-9 of "Cormen et al: Introduction to
 * agorithms", Chapter 22.5. 
 */
void generateDigraph1() {
	digraph = new Digraph();

	// creat Node
	Node* a = new Node("a", 1, 10);
	Node* b = new Node("b", 2, 10);
	Node* c = new Node("c", 3, 10);
	Node* d = new Node("d", 4, 10);
	Node* e = new Node("e", 5, 10);
	Node* f = new Node("f", 6, 10);
	Node* g = new Node("g", 7, 10);
	Node* h = new Node("h", 8, 10);

	digraph->add_node(a);
	digraph->add_node(b);
	digraph->add_node(c);
	digraph->add_node(d);
	digraph->add_node(e);
	digraph->add_node(f);
	digraph->add_node(g);
	digraph->add_node(h);

	// Create Edge
	Edge* ab = new Edge(a,b);
	Edge* bf = new Edge(b,f);
	Edge* be = new Edge(b,e);
	Edge* ea = new Edge(e,a);
	Edge* ef = new Edge(e,f);
	Edge* bc = new Edge(b,c);
	Edge* fg = new Edge(f,g);
	Edge* gf = new Edge(g,f);
	Edge* cd = new Edge(c,d);
	Edge* dc = new Edge(d,c);
	Edge* cg = new Edge(c,g);
	Edge* gh = new Edge(g,h);
	Edge* hh = new Edge(h,h);
	Edge* dh = new Edge(d,h);

	digraph->add_edge(ab);
	digraph->add_edge(bf);
	digraph->add_edge(be);
	digraph->add_edge(ea);
	digraph->add_edge(ef);
	digraph->add_edge(bc);
	digraph->add_edge(fg);
	digraph->add_edge(gf);
	digraph->add_edge(cd);
	digraph->add_edge(dc);
	digraph->add_edge(cg);
	digraph->add_edge(gh);
	digraph->add_edge(hh);
	digraph->add_edge(dh);

	if (output) {
		digraph->write_graphviz(std::cout);
	}
}

/**
 * An example of the paper RTNS2015.
 * Used to calculate the linear factor (or utilization)
 */
void generateDigraph2() {
	digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,2,3);
	Node* v2 = new Node("v2",1,2,3);
	Node* v3 = new Node("v3",1,1,2);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);

	// create edges
	Edge* v1_v1 = new Edge(v1,v1);
	v1_v1->set_separation_time(3);

	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(4);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(3);

	Edge* v3_v2 = new Edge(v3,v2);
	v3_v2->set_separation_time(2);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(2);

	digraph->add_edge(v1_v1);
	digraph->add_edge(v1_v2);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v3_v2);
	digraph->add_edge(v3_v1);

	if (output) {
		digraph->write_graphviz(std::cout);
	}
}

/**
 * An example of the paper RTS2014.
 * Used to test rbf and dbf
 */
void generateDigraph3() {
	digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,1,10);
	Node* v2 = new Node("v2",1,2,10);
	Node* v3 = new Node("v3",1,1,10);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);

	// create edges
	Edge* v1_v1 = new Edge(v1,v1);
	v1_v1->set_separation_time(10);

	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(20);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(10);

	Edge* v3_v2 = new Edge(v3,v2);
	v3_v2->set_separation_time(20);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(10);

	digraph->add_edge(v1_v1);
	digraph->add_edge(v1_v2);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v3_v2);
	digraph->add_edge(v3_v1);

	if (output) {
		digraph->write_graphviz(std::cout);
	}
}

/**
 * An example of Martin Stigge et al., The digraph real-time task model, RTAS2011
 * Used to dbf
 */
void generateDigraph4() {
	digraph = new Digraph();

	// create nodes
	Node* v1 = new Node("v1",1,2,5);
	Node* v2 = new Node("v2",1,1,8);
	Node* v3 = new Node("v3",1,3,8);
	Node* v4 = new Node("v4",1,5,10);
	Node* v5 = new Node("v5",1,1,5);

	digraph->add_node(v1);
	digraph->add_node(v2);
	digraph->add_node(v3);
	digraph->add_node(v4);
	digraph->add_node(v5);

	// create edges
	Edge* v1_v2 = new Edge(v1,v2);
	v1_v2->set_separation_time(10);

	Edge* v1_v5 = new Edge(v1,v5);
	v1_v5->set_separation_time(20);

	Edge* v2_v3 = new Edge(v2,v3);
	v2_v3->set_separation_time(15);

	Edge* v2_v4 = new Edge(v2,v4);
	v2_v4->set_separation_time(20);

	Edge* v3_v1 = new Edge(v3,v1);
	v3_v1->set_separation_time(11);

	Edge* v4_v2 = new Edge(v4,v2);
	v4_v2->set_separation_time(20);

	Edge* v5_v4 = new Edge(v5,v4);
	v5_v4->set_separation_time(10);

	digraph->add_edge(v1_v2);
	digraph->add_edge(v1_v5);
	digraph->add_edge(v2_v3);
	digraph->add_edge(v2_v4);
	digraph->add_edge(v3_v1);
	digraph->add_edge(v4_v2);
	digraph->add_edge(v5_v4);

	if (output) {
		digraph->write_graphviz(std::cout);
	}
}

// Compare two matrices
bool compare_tow_matrices(double A[][9], double** B, int n) {
	for (int i=0; i<n; i++) for (int j=0; j<n; j++)
		if (abs(A[i][j]-B[i][j])>EPSILON) return false;
	return true;
}

bool compare_tow_matrices(double** A, double** B, int n) {
	for (int i=0; i<n; i++) for (int j=0; j<n; j++)
		if (abs(A[i][j]-B[i][j])>EPSILON) {
			if (output) {
				cout<<i<<","<<j<<endl;
				cout<<A[i][j]<<"..."<<B[i][j]<<endl;
			}
			return false;
		}
	return true;
}

TEST(DigraphTest, DeepFirstSearch)
{
	generateDigraph0();

	digraph->generate_strongly_connected_components();

	vector<Digraph*> sccs = digraph->sccs;

	EXPECT_EQ(sccs.size(),4);

	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}
}

/**
 * We test the constructor function and algorithm of calculating strongly connected components of Digraph class.
 * The digraph comes from Figure 22-9 of "Cormen et al: Introduction to
 * agorithms", Chapter 22.5.
 */
TEST(DigraphTest, StronglyConnectedComponents)
{
	generateDigraph1();

	digraph->generate_strongly_connected_components();

	vector<Digraph*> sccs = digraph->sccs;

	EXPECT_EQ(sccs.size(),4);
	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}
}

/**
 * Test the generations of UDRT and GDRT
 */
TEST(DigraphTest, TransformedDigraph)
{
	generateDigraph2();

	digraph->calculate_gcd();

	digraph->prepare_rbf_calculation(false);
	digraph->prepare_ibf_calculation(false);

	if (output) {
		digraph->unit_digraph->write_graphviz(cout);
		digraph->gran_digraph->write_graphviz(cout);
	}
}

/**
 * Test karp's algorithm for calculating the maximum cycle mean
 */
TEST(DigraphTest, Digraph2)
{
	generateDigraph2();

	digraph->generate_strongly_connected_components();
	vector<Digraph*> sccs = digraph->sccs;
	EXPECT_EQ(sccs.size(),1);

	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}

	digraph->check_strongly_connected();
	EXPECT_EQ(digraph->strongly_connected,true);

	digraph->calculate_gcd();
	EXPECT_EQ(digraph->gcd,1);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,2.0/3);

	digraph->calculate_linear_upper_bounds();
	EXPECT_EQ(digraph->c_sum,5);
	EXPECT_EQ((int)(digraph->c_rbf*1000), 2000);
	EXPECT_EQ((int)(digraph->c_ibf*1000), 666);
	EXPECT_EQ((int)(digraph->c_dbf*1000), 666);

	digraph->tf = 100;

	digraph->prepare_rbf_calculation(false);

	double matrix[9][9] = {{0,2,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{0,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,2,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,1},
	{0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0}};

	bool test = compare_tow_matrices(matrix,digraph->unit_digraph->matrix,9);

	Utility::output_matrix(digraph->unit_digraph->matrix, digraph->unit_digraph->n_size, digraph->unit_digraph->n_size);
	EXPECT_EQ(test, true);
	EXPECT_EQ(digraph->unit_digraph->lper, 3);

	//cout<<"linear defect = " << digraph->unit_digraph->ldef<<endl;

	// test (ldef,ldef+10]
	map<int,double**> matrices;
	matrices[1] = digraph->unit_digraph->matrix;
	int ldef = MaxPlusAlgebra::calculate_linear_defect(matrices,9,2.0/3,3,100);

	EXPECT_EQ(digraph->unit_digraph->ldef, ldef);
	
	for (int t=2; t<=100; t++) {
		matrices[t] = MaxPlusAlgebra::multiply_maxplus_matrix(matrices[(t-1)],matrices[1],9);
	}

	for (int t=digraph->unit_digraph->ldef+1; t<=50; t++) {
		//cout<<"=================="<<t<<"====================="<<endl;

		double** temp = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(matrices[t],9,digraph->linear_factor,digraph->unit_digraph->lper);
		double** temp2 = matrices[t+digraph->unit_digraph->lper]; 
		
		EXPECT_EQ(compare_tow_matrices(temp,temp2,9), true);
	}

	// test rbf
	double rbf[16] = {0,2,2,3,4,4,5,6,6,7,8,8,9,10,10,11};
	for (int i=0; i<16;i++)
		EXPECT_EQ(digraph->rbf(i),rbf[i]);

	digraph->prepare_ibf_calculation(false);

	EXPECT_EQ(digraph->unit_digraph->lper,digraph->gran_digraph->lper);
	EXPECT_EQ(digraph->unit_digraph->ldef,digraph->gran_digraph->ldef);
	// test ibf

	// Output maximum element
	if (false) {
		for ( int i=0; i<20; i++)
			cout<<digraph->gran_digraph->maximum_element_map[i]<<endl;
	}

	double ibf[16] = {0,1,2,2,3,4,4,5,6,6,7,8,8,9,10,10};
	for (int i=0; i<16;i++)
		EXPECT_EQ(digraph->ibf(i),ibf[i]);

	// test dbf
	
	if (output) {
		cout<<digraph->unit_digraph->iSet.size()<<endl;
		for ( int i=0; i<20; i++)
			cout<<digraph->dbf(i)<<endl;
	}
}

TEST(DigraphTest, Digraph0)
{
	generateDigraph0();
	
	digraph->generate_strongly_connected_components();
    vector<Digraph*> sccs = digraph->sccs;
	EXPECT_EQ(sccs.size(),4);

	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}

	digraph->check_strongly_connected();
	EXPECT_EQ(digraph->strongly_connected,false);

	digraph->calculate_gcd();
	EXPECT_EQ(digraph->gcd,4);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,0.75);

}
	
/**
 * An example of the paper RTS2014.
 * Used to test rbf and dbf
 */
TEST(DigraphTest, Digraph3)
{
	generateDigraph3();

	digraph->generate_strongly_connected_components();
	vector<Digraph*> sccs = digraph->sccs;
	EXPECT_EQ(sccs.size(),1);

	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}

	digraph->check_strongly_connected();
	EXPECT_EQ(digraph->strongly_connected,true);

	digraph->calculate_gcd();
	EXPECT_EQ(digraph->gcd,10);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,0.1);

	digraph->tf = 100;

	digraph->prepare_rbf_calculation(false);

	// digraph->unit_digraph->write_graphviz(cout);

	EXPECT_EQ(digraph->unit_digraph->lper*digraph->gcd, 10);

	if(output) {
		Utility::output_matrix(digraph->unit_digraph->matrix,digraph->unit_digraph->n_size,digraph->unit_digraph->n_size);
	}

	EXPECT_EQ(digraph->unit_digraph->lper,1);
	EXPECT_EQ(digraph->unit_digraph->ldef,6);

	if(output) {
		digraph->unit_digraph->write_graphviz(cout);
	}

	//EXPECT_EQ(digraph->unit_digraph->ldef, 6);

	// test rbf
	if (output) {
		double rbf[50];
		for (int i=0; i<50; i++) rbf[i] = (double)i/10;

		for (int i=80; i<1000;i+=10) {
			cout<<digraph->rbf(i)<<endl;
			cout<<digraph->dbf(i)<<endl;
		}
	}

	double** A = new double*[5];
	for (int i=0; i<5; i++) A[i] = new double[5];
	for (int i=0; i<5; i++) for (int j=0; j<5; j++) A[i][j] = NEG_INFINITY;

	A[0][0] = 1;
	A[0][3] = 1;
	A[1][1] = 0;
	A[1][2] = 2;
	A[2][0] = 1;
	A[2][2] = 0;
	A[2][4] = 1;
	A[3][1] = 0;
	A[4][1] = 0;

	std::map<int,double**> matrices;
	matrices[1]=A;

	int ldef = MaxPlusAlgebra::calculate_linear_defect(matrices, 5, 1, 1, 100);
	//cout<<digraph->unit_digraph->lfac<<","<<digraph->unit_digraph->lper<<endl;
	EXPECT_EQ(ldef,6);

	// test ibf

	digraph->prepare_ibf_calculation(false);

	EXPECT_EQ(digraph->unit_digraph->lper*digraph->unit_digraph->gcd,digraph->gran_digraph->lper*digraph->gran_digraph->gcd);
	//EXPECT_EQ(digraph->unit_digraph->ldef,digraph->gran_digraph->ldef);

	

	if (output) {
		digraph->gran_digraph->write_graphviz(cout);
		for (int i=0; i<20;i++) {
			cout<<"=================="<<i<<"============="<<endl;
			//Utility::output_matrix(digraph->unit_digraph->matrix_map[8],digraph->unit_digraph->n_size,digraph->unit_digraph->n_size);
			//cout<<digraph->unit_digraph->maximum_element_map[8]<<endl;
			cout<<digraph->rbf(i)<<endl;
			cout<<digraph->ibf(i)<<endl;
			cout<<digraph->dbf(i)<<endl;
			
		}
	}

	int rbf[15] = {0,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	for (int i=0; i<15; i++) 
		EXPECT_EQ(digraph->rbf(i*digraph->gcd),rbf[i]);
	
}

 TEST(DigraphTest, Digraph4)
 {
	 generateDigraph4();
	 digraph->tf = 100;

	 digraph->prepare_digraph();
	 
	 digraph->prepare_rbf_calculation(false);

	 if (output) {
		for (int i=0; i<60;i++) {
			//cout<<digraph->rbf(i)<<endl;
			cout<<"dbf("<<i<<")="<<digraph->dbf(i)<<endl;
		}
	}
 }

 /**
 * Test karp's algorithm for calculating the maximum cycle mean
 */
TEST(DigraphTest, Digraph5)
{
	generateDigraph2();

	for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
		Node* node = *iter;
		node->wcet = node->wcet*1000;
		node->deadline = node->deadline*1000;
	}

	for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
		Edge* edge = *iter;
		edge->separationTime = edge->separationTime*1000;
	}


	digraph->generate_strongly_connected_components();
	vector<Digraph*> sccs = digraph->sccs;
	EXPECT_EQ(sccs.size(),1);

	if (output) {
		for (vector<Digraph*>::iterator iter = sccs.begin(); iter != sccs.end(); iter++) {
			(*iter)->write_graphviz(std::cout);
		}
	}

	digraph->check_strongly_connected();
	EXPECT_EQ(digraph->strongly_connected,true);

	digraph->calculate_gcd();
	EXPECT_EQ(digraph->gcd,1000);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,2.0/3);

	digraph->calculate_linear_upper_bounds();
	EXPECT_EQ(digraph->c_sum,5000);
	EXPECT_EQ((int)(digraph->c_rbf), 2000);
	EXPECT_EQ((int)(digraph->c_ibf), 666);
	EXPECT_EQ((int)(digraph->c_dbf), 666);

	digraph->tf = 100;

	digraph->prepare_rbf_calculation(false);

	double matrix[9][9] = {{0,2,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{0,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,2,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,0,NEG_INFINITY},
	{NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,1},
	{0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0,NEG_INFINITY,NEG_INFINITY,NEG_INFINITY,0}};

	bool test = compare_tow_matrices(matrix,digraph->unit_digraph->matrix,9);

	//Utility::output_matrix(digraph->unit_digraph->matrix, digraph->unit_digraph->n_size, digraph->unit_digraph->n_size);
	//EXPECT_EQ(test, true);
	EXPECT_EQ(digraph->unit_digraph->lper, 3);

	//cout<<"linear defect = " << digraph->unit_digraph->ldef<<endl;

	// test (ldef,ldef+10]
	map<int,double**> matrices;
	matrices[1] = digraph->unit_digraph->matrix;
	int ldef = MaxPlusAlgebra::calculate_linear_defect(matrices,9,2.0/3,3,100);

	EXPECT_EQ(digraph->unit_digraph->ldef, ldef);
}
#endif