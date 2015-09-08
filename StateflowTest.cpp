#include <gtest/gtest.h>

#include "Stateflow.h"

Stateflow* sf;

bool output = false;
extern double EPSILON;

bool compare_two_matrices(double** A, double** B, int n) {
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

/**
 * An example for Figure 1 in 
 * Zeng and Di Natale, Schedulability analysis of periodic tasks implementing synchronous finite state machines.
 * ECRTS2012.
 */
void generateStateflow0() {
	sf = new Stateflow(100);

	State* s1 = new State("s1",1);
	State* s2 = new State("s2",2);
	State* s3 = new State("s3",3);

	sf->add_state(s1);
	sf->add_state(s2);
	sf->add_state(s3);

	Transition* t12 = new Transition(s1,s2);
	t12->wcet = 25;
	t12->period = 200;
	t12->priority = 1;
	sf->add_transition(t12);

	Transition* t23_0 = new Transition(s2,s3);
	t23_0->wcet = 10;
	t23_0->period = 200;
	t23_0->priority = 1;
	sf->add_transition(t23_0);

	Transition* t23_1 = new Transition(s2,s3);
	t23_1->wcet = 15;
	t23_1->period = 500;
	t23_1->priority = 2;
	sf->add_transition(t23_1);

	Transition* t31 = new Transition(s3,s1);
	t31->wcet = 30;
	t31->period = 500;
	t31->priority = 2;
	sf->add_transition(t31);

	if (output) sf->write_graphviz(cout);
}

TEST(StateflowTest, Stateflow0)
{
	// rbf_{i,j}[s,f)
	double**** rbfijsf = new double***[3];
	for (int i=0; i<3; i++) rbfijsf[i] = new double**[3];
	for (int i=0; i<3; i++) for (int j=0; j<3; j++) rbfijsf[i][j] = new double*[7];
	for (int i=0; i<3; i++) for (int j=0; j<3; j++) for (int s=0; s<7; s++) rbfijsf[i][j][s] = new double[7];
	for (int i=0; i<3; i++) for (int j=0; j<3; j++) 
		for (int s=0; s<7; s++) for (int f=0; f<7; f++)
			rbfijsf[i][j][s][f] = NEG_INFINITY;

	// s1->s1
	rbfijsf[0][0][0][3] = rbfijsf[0][0][0][2] = rbfijsf[0][0][0][1] = rbfijsf[0][0][0][0] = 0;
	rbfijsf[0][0][0][6] = rbfijsf[0][0][0][5] = rbfijsf[0][0][0][4] = 65;

	rbfijsf[0][0][1][3] = rbfijsf[0][0][1][2] = rbfijsf[0][0][1][1] = 0;
	rbfijsf[0][0][1][6] = rbfijsf[0][0][1][5] = rbfijsf[0][0][1][4] = 65;

	rbfijsf[0][0][2][6] = rbfijsf[0][0][2][5] = rbfijsf[0][0][2][4] = rbfijsf[0][0][2][3] = rbfijsf[0][0][2][2] = 0;

	rbfijsf[0][0][3][6] = rbfijsf[0][0][3][5] = rbfijsf[0][0][3][4] = rbfijsf[0][0][3][3] = 0;

	rbfijsf[0][0][4][6] = rbfijsf[0][0][4][5] = rbfijsf[0][0][4][4] = 0;

	rbfijsf[0][0][5][6] = rbfijsf[0][0][5][5] = 0;
	
	rbfijsf[0][0][6][6] = 0;

	// s1->s2
	rbfijsf[0][1][0][4] = rbfijsf[0][1][0][3] = rbfijsf[0][1][0][2] = rbfijsf[0][1][0][1] = 25;
	rbfijsf[0][1][0][6] = rbfijsf[0][1][0][5] = 90;

	rbfijsf[0][1][1][4] = rbfijsf[0][1][1][3] = rbfijsf[0][1][1][2] = 25;
	rbfijsf[0][1][1][6] = rbfijsf[0][1][1][5] = 90;

	rbfijsf[0][1][2][6] = rbfijsf[0][1][2][5] = rbfijsf[0][1][2][4] = rbfijsf[0][1][2][3] = 25;

	rbfijsf[0][1][3][6] = rbfijsf[0][1][3][5] = 25;

	rbfijsf[0][1][4][6] = rbfijsf[0][1][4][5] = 25;

	rbfijsf[0][1][5][6] = 25;

	// s1->s3
	rbfijsf[0][2][0][3] = rbfijsf[0][2][0][2] = 35;
	rbfijsf[0][2][0][5] = rbfijsf[0][2][0][4] = 40;
	rbfijsf[0][2][0][6] = 100;

	rbfijsf[0][2][1][3] = 35;
	rbfijsf[0][2][1][5] = rbfijsf[0][2][1][4] = 40;
	rbfijsf[0][2][1][6] = 100;

	rbfijsf[0][2][2][6] = rbfijsf[0][2][2][5] = rbfijsf[0][2][2][4] = 40;

	rbfijsf[0][2][3][6] = 35;

	rbfijsf[0][2][4][6] = 35;

	// s2->s1
	rbfijsf[1][0][0][6] = rbfijsf[1][0][0][5] = rbfijsf[1][0][0][4] = 45;

	rbfijsf[1][0][1][6] = rbfijsf[1][0][1][5] = rbfijsf[1][0][1][4] = 40;

	rbfijsf[1][0][2][6] = rbfijsf[1][0][2][5] = rbfijsf[1][0][2][4] = 40;

	// s2->s2
	rbfijsf[1][1][0][4] = rbfijsf[1][1][0][3] = rbfijsf[1][1][0][2] = rbfijsf[1][1][0][1] = rbfijsf[1][1][0][0] = 0;
	rbfijsf[1][1][0][6] = rbfijsf[1][1][0][5] = 70;

	rbfijsf[1][1][1][4] = rbfijsf[1][1][1][3] = rbfijsf[1][1][1][2] = rbfijsf[1][1][1][1] = 0;
	rbfijsf[1][1][1][6] = rbfijsf[1][1][1][5] = 65;

	rbfijsf[1][1][2][4] = rbfijsf[1][1][2][3] = rbfijsf[1][1][2][2] = 0;
	rbfijsf[1][1][2][6] = rbfijsf[1][1][2][5] = 65; 

	rbfijsf[1][1][3][6] = rbfijsf[1][1][3][5] = rbfijsf[1][1][3][4] = rbfijsf[1][1][3][3] = 0;

	rbfijsf[1][1][4][6] = rbfijsf[1][1][4][5] = rbfijsf[1][1][4][4] = 0;

	rbfijsf[1][1][5][6] = rbfijsf[1][1][5][5] = 0;
	
	rbfijsf[1][1][6][6] = 0;

	// s2->s3
	rbfijsf[1][2][0][5] = rbfijsf[1][2][0][4] = rbfijsf[1][2][0][3] = rbfijsf[1][2][0][2] = rbfijsf[1][2][0][1] = 15;
	rbfijsf[1][2][0][6] = 80;

	rbfijsf[1][2][1][3] = rbfijsf[1][2][1][2] = 10;
	rbfijsf[1][2][1][5] = rbfijsf[1][2][1][4] = 15;
	rbfijsf[1][2][1][6] = 75;

	rbfijsf[1][2][2][3] = 10;
	rbfijsf[1][2][2][5] = rbfijsf[1][2][2][4] = 15;
	rbfijsf[1][2][2][6] = 75;

	rbfijsf[1][2][3][6] = rbfijsf[1][2][3][5] = rbfijsf[1][2][3][4] = 15;

	rbfijsf[1][2][4][6] = rbfijsf[1][2][4][5] = 10;

	rbfijsf[1][2][5][6] = 10;

	// s3->s1
	rbfijsf[2][0][0][3] = rbfijsf[2][0][0][2] = rbfijsf[2][0][0][1] = 30;
	rbfijsf[2][0][0][6] = rbfijsf[2][0][0][5] = rbfijsf[2][0][0][4] = 95;

	rbfijsf[2][0][1][6] = rbfijsf[2][0][1][5] = rbfijsf[2][0][1][4] = 30;

	rbfijsf[2][0][2][6] = rbfijsf[2][0][2][5] = rbfijsf[2][0][2][4] = 30;

	rbfijsf[2][0][3][6] = rbfijsf[2][0][3][5] = rbfijsf[2][0][3][4] = 30;

	// s3->s2
	rbfijsf[2][1][0][4] = rbfijsf[2][1][0][3] = rbfijsf[2][1][0][2] = 55;
	rbfijsf[2][1][0][6] = rbfijsf[2][1][0][5] = 120;

	rbfijsf[2][1][1][6] = rbfijsf[2][1][1][5] = 55;

	rbfijsf[2][1][2][6] = rbfijsf[2][1][2][5] = 55;

	rbfijsf[2][1][3][6] = rbfijsf[2][1][3][5] = 55;

	// s3->s3
	rbfijsf[2][2][0][2] = rbfijsf[2][2][0][1] = rbfijsf[2][2][0][0] = 0;
	rbfijsf[2][2][0][3] = 65;
	rbfijsf[2][2][0][5] = rbfijsf[2][2][0][4] = 70;
	rbfijsf[2][2][0][6] = 130;

	rbfijsf[2][2][1][5] = rbfijsf[2][2][1][4] = rbfijsf[2][2][1][3] = rbfijsf[2][2][1][2] = rbfijsf[2][2][1][1] = 0;
	rbfijsf[2][2][1][6] = 65;

	rbfijsf[2][2][2][5] = rbfijsf[2][2][2][4] = rbfijsf[2][2][2][3] = rbfijsf[2][2][2][2] = 0;
	rbfijsf[2][2][2][6] = 65;

	rbfijsf[2][2][3][5] = rbfijsf[2][2][3][4] = rbfijsf[2][2][3][3] = 0;
	rbfijsf[2][2][3][6] = 65;

	rbfijsf[2][2][4][6] = rbfijsf[2][2][4][5] = rbfijsf[2][2][4][4] = 0;

	rbfijsf[2][2][5][6] = rbfijsf[2][2][5][5] = 0;

	rbfijsf[2][2][6][6] = 0;

	// Execution Request Matrix
	double** exec_matrix = new double*[3];
	for (int i=0; i<3; i++) exec_matrix[i] = new double[3];
	exec_matrix[0][0] = 65;
	exec_matrix[0][1] = 90;
	exec_matrix[0][2] = 100;
	exec_matrix[1][0] = 45;
	exec_matrix[1][1] = 70;
	exec_matrix[1][2] = 80;
	exec_matrix[2][0] = 95;
	exec_matrix[2][1] = 120;
	exec_matrix[2][2] = 130;


	generateStateflow0();

	sf->calculate_gcd();
	EXPECT_EQ(sf->gcd,100);

	sf->calculate_t_gcd();
	EXPECT_EQ(sf->t_gcd,5);

	sf->calculate_hyperperiod();
	EXPECT_EQ(sf->hyperperiod,1000);

	sf->set_state_number();
	EXPECT_EQ(sf->n_state,3);

	sf->generate_rbf_time_instances();
	if (output) {
		for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter++) {
			cout<<*iter<<"\t";
		}
		cout<<endl;

		cout<<sf->index_time.size()<<endl;
		cout<<sf->time_index.size()<<endl;
	}

	sf->generate_rbfs();
	for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
		EXPECT_EQ(compare_two_matrices(sf->rbfs[i][j], rbfijsf[i][j], 7), true);
	}

	if (output) {
		for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
			cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
			Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
		}
	}

	sf->generate_exec_req_matrix();
	EXPECT_EQ(compare_two_matrices(sf->exec_req_matrix,exec_matrix, 3), true);
	
	if (output) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);

	sf->generate_simple_digraph();
	if (output) {
		sf->simple_digraph->write_graphviz(cout);
		cout<<sf->simple_digraph->linear_factor<<"\t"<<sf->simple_digraph->c_sum<<"\t"<<sf->simple_digraph->c_rbf<<"\t"<<sf->simple_digraph->c_ibf
			<<"\t"<<sf->simple_digraph->c_dbf<<endl;
	}

	sf->generate_precise_digraph();
	if (output) {
		sf->precise_digraph->write_graphviz(cout);
		cout<<sf->precise_digraph->linear_factor<<"\t"<<sf->precise_digraph->c_sum<<"\t"<<sf->precise_digraph->c_rbf<<"\t"<<sf->precise_digraph->c_ibf
			<<"\t"<<sf->precise_digraph->c_dbf<<endl;
	}

	sf->calculate_linear_factor();

	EXPECT_EQ(sf->lfac, sf->precise_digraph->linear_factor);

	sf->check_irreducible();
	EXPECT_EQ(sf->isIrred, true);

	sf->calculate_generialized_period();
	EXPECT_EQ(sf->gper,1);

	sf->tf0 = 100;

	sf->calculate_generialized_defect();
	EXPECT_EQ(sf->gdef,1);

	sf->calculate_exec_req_matrix_power(10);

	if (output) {
		for (int i=2; i<10; i++) {
			cout<<"================="<<i<<"================="<<endl;
			Utility::output_matrix(sf->exec_req_matrix_power[i],3,3);
		}
	}

	// test rbf[s,f)
	EXPECT_EQ(sf->get_rbf(0,800), 120);
	EXPECT_EQ(sf->get_rbf(0,900), 130);
	EXPECT_EQ(sf->get_rbf(0,1000), 130);

	EXPECT_EQ(sf->get_rbf(200,800), 90);
	EXPECT_EQ(sf->get_rbf(400,900), 75);
	EXPECT_EQ(sf->get_rbf(300,1000), 75);

	EXPECT_EQ(sf->get_rbf(1200,1800), 90);
	EXPECT_EQ(sf->get_rbf(2400,2900), 75);
	EXPECT_EQ(sf->get_rbf(3300,4000), 75);

	EXPECT_EQ(sf->get_rbf(200,1200), 130);
	EXPECT_EQ(sf->get_rbf(400,1500), 140);

	if (output) {
		for (int i=0; i<2000; i+=100) {
			for (int j=0; j<2000; j+=100)
				cout<<"rbf["<<+i<<","<<j<<")="<<sf->get_rbf(i,j)<<endl;
		}
	}

	// test rbf(t)
	sf->simple_digraph->tf = 2000;

	sf->simple_digraph->prepare_rbf_calculation(false);
	if (output) {
		for (int i=0; i<2000; i+=100)
			cout<<"rbf("<<i<<")="<<sf->simple_digraph->rbf(i)<<endl;
	}

	sf->precise_digraph->tf = 1000;

	sf->precise_digraph->prepare_rbf_calculation(false);

	if (output) {
		int s = 0;
		for (int f = 0; f<1000; f+=100)
			cout<<"rbf["<<+s<<","<<f<<")="<<sf->get_rbf(s,f)<<endl;

		for (int t = 100; t<1000; t+=100) {
			cout<<"Stateflow:rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
			cout<<"Digraph:rbf("<<t<<")="<<sf->precise_digraph->rbf(t)<<endl;
		}

		cout<<"Stateflow:rbf("<<900<<")="<<sf->get_rbf(900)<<endl;
	}

	for (int t=0; t<10000; t+=100)
		EXPECT_EQ(sf->get_rbf(t),sf->precise_digraph->rbf(t));

	// test dbf[s,f)
	if (output) {
		for (int i=0; i<100; i+=100) {
			for (int j=0; j<2000; j+=100) {
				cout<<"rbf["<<+i<<","<<j<<")="<<sf->get_rbf(i,j)<<endl;
				cout<<"dbf["<<+i<<","<<j<<")="<<sf->get_dbf(i,j)<<endl;
			}
		}
	}

	for (int s=0; s<10000; s+=100) 	for (int f=0; f<10000; f+=100) {
		if (sf->get_rbf(s,f)!=sf->get_ibf(s,f))
			cout<<"Error goes here:"<<s<<","<<f<<endl;
		EXPECT_EQ(sf->get_rbf(s,f),sf->get_ibf(s,f));	
		EXPECT_EQ(sf->get_rbf(s,f),sf->get_dbf(s,f));
	}

	// test dbf(t)
	for (int t=0; t<10000; t+=100)
		EXPECT_EQ(sf->get_rbf(t),sf->get_dbf(t));

	// test ibf[s,f)
	if (output) {
		for (int i=0; i<100; i+=100) {
			for (int j=0; j<1000; j+=10) {
				cout<<"rbf["<<+i<<","<<j<<")="<<sf->get_rbf(i,j)<<endl;
				cout<<"ibf["<<+i<<","<<j<<")="<<sf->get_ibf(i,j)<<endl;
				cout<<"dbf["<<+i<<","<<j<<")="<<sf->get_dbf(i,j)<<endl;
			}
		}
	}

	for (int i=0; i<1000; i+=10) for (int j=0; j<1000; j+=10) 
		EXPECT_LE(sf->get_ibf(i,j),sf->get_rbf(i,j));

	// test ibf(t)
	for (int t=0; t<10000; t+=100)
		EXPECT_EQ(sf->get_rbf(t),sf->get_ibf(t));

	for (int t=0; t<1000; t+=10)
		EXPECT_GE(sf->get_rbf(t),sf->get_ibf(t));

	if (output) {
		for (int t=0; t<1000; t+=10) {
			cout<<"ibf("<<+t<<")="<<sf->get_ibf(t)<<endl;
		}
	}

}