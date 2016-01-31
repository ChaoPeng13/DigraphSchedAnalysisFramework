#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "Digraph.h"
#include "MaxPlusAlgebra.h"
#include "DigraphExample.h"

Digraph* digraph;

#define output false
extern double EPSILON;

TEST(DigraphTest, DeepFirstSearch)
{
	Digraph* digraph = DigraphExample::generateDigraph0();

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
	Digraph* digraph = DigraphExample::generateDigraph1();

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
	Digraph* digraph = DigraphExample::generateDigraph2();

	digraph->calculate_period_gcd();
	digraph->calculate_all_gcd();

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
	Digraph* digraph = DigraphExample::generateDigraph2();

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

	digraph->calculate_period_gcd();
	digraph->calculate_all_gcd();
	EXPECT_EQ(digraph->pGCD,1);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,2.0/3);

	digraph->calculate_csum();
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

	double** matrix2 = new double*[9];
	for (int i=0; i<9; i++) matrix2[i] = new double[9];
	for (int i=0; i<9; i++) for (int j=0; j<9; j++) matrix2[i][j] = matrix[i][j];

	bool test = Utility::compare_two_matrices(matrix2,digraph->unit_digraph->matrix,9);

	if (output) Utility::output_matrix(digraph->unit_digraph->matrix, digraph->unit_digraph->n_size, digraph->unit_digraph->n_size);
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
		
		EXPECT_EQ(Utility::compare_two_matrices(temp,temp2,9), true);
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
	
	if (false) {
		cout<<digraph->unit_digraph->iSet.size()<<endl;
		for ( int i=0; i<20; i++)
			cout<<digraph->dbf(i)<<endl;
	}
}

TEST(DigraphTest, Digraph0)
{
	Digraph* digraph = DigraphExample::generateDigraph0();
	
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

	digraph->calculate_period_gcd();
	EXPECT_EQ(digraph->pGCD,4);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,0.75);

}
	
/**
 * An example of the paper RTS2014.
 * Used to test rbf and dbf
 */
TEST(DigraphTest, Digraph3)
{
	Digraph* digraph = DigraphExample::generateDigraph3();

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

	digraph->calculate_period_gcd();
	digraph->calculate_all_gcd();
	EXPECT_EQ(digraph->pGCD,10);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,0.1);

	digraph->tf = 1000;

	digraph->prepare_rbf_calculation(false);

	// digraph->unit_digraph->write_graphviz(cout);

	EXPECT_EQ(digraph->unit_digraph->lper*digraph->pGCD, 10);

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
	if (false) {
		double rbf[50];
		for (int i=0; i<50; i++) rbf[i] = (double)i/10;

		for (int i=0; i<1000;i+=10) {
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
		EXPECT_EQ(digraph->rbf(i*digraph->pGCD),rbf[i]);
	
}

 TEST(DigraphTest, Digraph4)
 {
	 Digraph* digraph = DigraphExample::generateDigraph4();
	 digraph->tf = 100;

	 digraph->prepare_digraph();
	 
	 digraph->prepare_rbf_calculation(false);

	 int dbf[60] = {0, 0, 0, 0, 0, 2, 2, 2, 3, 3,
					5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
					6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
					6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
					8, 8, 8, 9, 9, 9, 9, 9, 9, 9,
					11, 11, 11, 11, 11, 11, 11, 11, 11, 11};

	
	for (int i=0; i<60;i++) {
		EXPECT_EQ(digraph->dbf(i), dbf[i]);
	}
 }

 /**
 * Test karp's algorithm for calculating the maximum cycle mean
 */
TEST(DigraphTest, Digraph5)
{
	Digraph* digraph = DigraphExample::generateDigraph2();

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

	digraph->calculate_period_gcd();
	EXPECT_EQ(digraph->pGCD,1000);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,2.0/3);

	digraph->calculate_csum();
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

	double** matrix2 = new double*[9];
	for (int i=0; i<9; i++) matrix2[i] = new double[9];
	for (int i=0; i<9; i++) for (int j=0; j<9; j++) matrix2[i][j] = matrix[i][j];

	bool test = Utility::compare_two_matrices(matrix2,digraph->unit_digraph->matrix,9);

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