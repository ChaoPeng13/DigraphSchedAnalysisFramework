#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "Stateflow.h"
#include "StateflowExample.h"
#include "DigraphExample.h"
#include "MaxPlusAlgebra.h"

Stateflow* sf;

bool output = false;

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


	Stateflow* sf = StateflowExample::generateStateflow0();

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
		EXPECT_EQ(Utility::compare_two_matrices(sf->rbfs[i][j], rbfijsf[i][j], 7), true);
	}

	if (output) {
		for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
			cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
			Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
		}
	}

	sf->generate_exec_req_matrix();
	EXPECT_EQ(Utility::compare_two_matrices(sf->exec_req_matrix,exec_matrix, 3), true);
	
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
	sf->calculate_linear_upper_bounds(true);

	cout << "Linear factor = " << sf->lfac <<endl;
	cout << "cibf = " << sf->cibf << "\t crbf = " << sf->crbf <<endl;

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

	if (true) {
		for (int t=0; t<2000; t+=100) {
			cout << "rbf[0,"<<t<<")=" << sf->get_rbf(0,t) <<endl;
		}
	}

	// test rbf(t)
	sf->simple_digraph->tf = 2000;

	sf->simple_digraph->prepare_rbf_calculation(false);
	if (output) {
		for (int i=0; i<2000; i+=100)
			cout<<"rbf("<<i<<")="<<sf->simple_digraph->rbf(i)<<endl;
	}

	sf->precise_digraph->tf = 10000;

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

	// output all the rbfs for the example in Figure 1
	if (true) {
		int s = 0;
		cout << "rbf[s,f)"<<endl;
		for (int f = 0; f<1200; f+=100)
			cout<<"rbf["<<+s<<","<<f<<")="<<sf->get_rbf(s,f)<<endl;

		cout << "Simple Digraph rbf(t)"<<endl;
		for (int t = 0; t<1200; t+=100) {
			//cout<<"Stateflow:rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
			cout<<"Digraph:rbf("<<t<<")="<<sf->simple_digraph->rbf(t)<<endl;
		}

		cout << "Precise Digraph rbf(t)" <<endl;
		for (int t = 0; t<1200; t+=100) {
			cout<<"Stateflow:rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
			cout<<"Digraph:rbf("<<t<<")="<<sf->precise_digraph->rbf(t)<<endl;
		}

		cout<<"Stateflow:rbf("<<900<<")="<<sf->get_rbf(900)<<endl;
	}

	for (int t=0; t<10000; t+=100) {
		//cout<<t<<endl;
		EXPECT_EQ(sf->get_rbf(t),sf->precise_digraph->rbf(t));
	}

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

	if (true) {
		for (int t=0; t<=1000; t+=5) {
			cout << "ibf[0," << t << ")=" << sf->get_ibf(0,t) <<endl;
		}
	}

	if (false) {
		for (int t=0; t<=1000; t+=5) {
			cout << "ibf(" << t << ")=" << sf->get_ibf(t) <<endl;
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

	int nSize;
	cout << "lfac = " << sf->simple_digraph->linear_factor << endl; 
	cout << "lper = " << sf->simple_digraph->unit_digraph->lper << "\tldef = " << sf->simple_digraph->unit_digraph->ldef <<endl;

	/*
	double*** A = new double**[31];
	
	for (int t=1; t<=30; t++) {
		A[t] = sf->simple_digraph->rbf_exec_req_matrix(t,nSize);
	}

	for (int t=1; t<=26; t++) {
		bool equal = true;
		cout << "===========" << t << "===========" <<endl;
		Utility::output_matrix(A[t],nSize,nSize);
		cout << "+++++++++++" << t+4 << "+++++++++++" <<endl;
		Utility::output_matrix(A[t+4],nSize,nSize);
		for ( int i=0; i<nSize; i++) {
			for ( int j=0; j<nSize; j++) {
				if (abs(A[t][i][j] +65 - A[t+4][i][j])>0.000001) {
					equal = false;
					break;
				}
			}
			if (!equal) break;
		}
		cout << "t=" << t << "=>" << equal <<endl;
	}

	delete[] A;
	*/
	cout << "nSize of precise digraph's UDRT = " << sf->precise_digraph->unit_digraph->n_size <<endl;
	cout << "lfac = " << sf->precise_digraph->linear_factor << endl; 
	cout << "lper = " << sf->precise_digraph->unit_digraph->lper << "\tldef = " << sf->precise_digraph->unit_digraph->ldef <<endl;
	
	
	//map<int,double**> matrices;
	//matrices[1] = sf->precise_digraph->rbf_exec_req_matrix(1,nSize);
	int period = 1000;
	for (int t=1; t<=100; t++) {
		bool equal = true;
		int x = t*100;
		/*
		cout << "===========" << t << "===========" <<endl;
		Utility::output_matrix(B[t],nSize,nSize);
		cout << "+++++++++++" << t+10 << "+++++++++++" <<endl;
		Utility::output_matrix(B[t+10],nSize,nSize);
		*/
		//double** B1 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices,t,nSize);
		//double** B2 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices,t+10,nSize);

		double B1 = sf->precise_digraph->rbf(x);
		double B2 = sf->precise_digraph->rbf(x+period);
		cout << B1 <<"\t" << B2 <<endl;

		/*
		if (t==1000) {
			cout << "===========" << t << "===========" <<endl;
			Utility::output_matrix(B1,nSize,nSize);
			cout << "+++++++++++" << t+10 << "+++++++++++" <<endl;
			Utility::output_matrix(B2,nSize,nSize);
		}
		*/

		/*
		for ( int i=0; i<nSize; i++) {
			for ( int j=0; j<nSize; j++) {
				if (abs(B1[i][j] +130 - B2[i][j])>0.000001) {
					equal = false;
					break;
				}
			}
			if (!equal) break;
		}
		*/

		if (abs(B1+0.13*period - B2) > 0.000001) equal = false;
		cout << "t=" << t << "=>" << equal <<endl;
	}
}

TEST(StateflowTest, Stateflow3)
{
	Stateflow* sf = StateflowExample::generateStateflow3();
	sf->calculate_gcd();
	sf->calculate_t_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->tf0 = 100;

	sf->generate_exec_req_matrix();
	if (true) {
		cout << "Execution Request Matrix for the FSM:" <<endl;
		Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	}

	sf->generate_simple_digraph();
	sf->generate_precise_digraph();

	sf->calculate_linear_factor();
	sf->calculate_linear_upper_bounds(true);

	cout << "Linear factor = " << sf->lfac <<endl;
	cout << "cibf = " << sf->cibf << "\t crbf = " << sf->crbf <<endl;

	sf->check_irreducible();
	sf->calculate_generialized_period();
	sf->calculate_generialized_defect();

	if (true) {
		cout << "gper = " << sf->gper << "\t gdef = " << sf->gdef <<endl;
	}

	sf->calculate_exec_req_matrix_power(10);

	if (true) {
		for (int t=0; t<=1000; t+=100) {
			cout << "rbf[0,"<<t<<")=" << sf->get_rbf(0,t) <<endl;
		}
	}

	sf->simple_digraph->tf = 2000;
	sf->simple_digraph->prepare_rbf_calculation(false);

	sf->precise_digraph->tf = 10000;
	sf->precise_digraph->prepare_rbf_calculation(false);

	// output execution request matrix of the inaccurate digraph model
	if (true) {
		cout << "Execution request matrix of the inaccurate digraph model:" << endl;
		Utility::output_matrix(sf->simple_digraph->unit_digraph->matrix,sf->simple_digraph->unit_digraph->n_size,sf->simple_digraph->unit_digraph->n_size);
	}

	// output all the rbfs for the example in Figure 1
	if (true) {
		int s = 0;
		cout << "rbf[s,f)"<<endl;
		for (int f = 0; f<1200; f+=100)
			cout<<"rbf["<<+s<<","<<f<<")="<<sf->get_rbf(s,f)<<endl;

		cout << "Simple Digraph rbf(t)"<<endl;
		for (int t = 0; t<1200; t+=100) {
			//cout<<"Stateflow:rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
			cout<<"Digraph:rbf("<<t<<")="<<sf->simple_digraph->rbf(t)<<endl;
		}

		cout << "Precise Digraph rbf(t)" <<endl;
		for (int t = 0; t<1200; t+=100) {
			cout<<"Stateflow:rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
			cout<<"Digraph:rbf("<<t<<")="<<sf->precise_digraph->rbf(t)<<endl;
		}

		for (int t = 0; t<=500; t+=5) {
			cout<<"Stateflow:ibf("<<t<<")="<<sf->get_ibf(t)<<endl;
			//cout<<"Digraph:ibf("<<t<<")="<<sf->precise_digraph->ibf(t)<<endl;
		}

		cout<<"Stateflow:rbf("<<900<<")="<<sf->get_rbf(900)<<endl;
	}

	if (true) {
		for (int t=200; t<=700; t+=5) {
			cout << "ibf[200," << t << ")=" << sf->get_ibf(200,t) <<endl;
			cout << "rbf[200," << t << ")=" << sf->get_rbf(200,t) <<endl;
		}
	}

	if (true) {
		for (int t=0; t<=500; t+=5) {
			cout << "ibf[0," << t << ")=" << sf->get_ibf(0,t) <<endl;
			cout << "rbf[0," << t << ")=" << sf->get_rbf(0,t) <<endl;
		}
	}

	int nSize;
	cout << "lfac = " << sf->simple_digraph->linear_factor << endl; 
	cout << "lper = " << sf->simple_digraph->unit_digraph->lper << "\tldef = " << sf->simple_digraph->unit_digraph->ldef <<endl;

	/*
	double*** A = new double**[31];
	
	for (int t=1; t<=30; t++) {
		A[t] = sf->simple_digraph->rbf_exec_req_matrix(t,nSize);
	}

	for (int t=1; t<=26; t++) {
		bool equal = true;
		cout << "===========" << t << "===========" <<endl;
		Utility::output_matrix(A[t],nSize,nSize);
		cout << "+++++++++++" << t+4 << "+++++++++++" <<endl;
		Utility::output_matrix(A[t+4],nSize,nSize);
		for ( int i=0; i<nSize; i++) {
			for ( int j=0; j<nSize; j++) {
				if (abs(A[t][i][j] +65 - A[t+4][i][j])>0.000001) {
					equal = false;
					break;
				}
			}
			if (!equal) break;
		}
		cout << "t=" << t << "=>" << equal <<endl;
	}

	delete[] A;
	*/
	cout << "nSize of precise digraph's UDRT = " << sf->precise_digraph->unit_digraph->n_size <<endl;
	cout << "lfac = " << sf->precise_digraph->linear_factor << endl; 
	cout << "lper = " << sf->precise_digraph->unit_digraph->lper << "\tldef = " << sf->precise_digraph->unit_digraph->ldef <<endl;
	
	
	//map<int,double**> matrices;
	//matrices[1] = sf->precise_digraph->rbf_exec_req_matrix(1,nSize);
	int period = 1000;
	for (int t=1; t<=100; t++) {
		bool equal = true;
		int x = t*100;
		/*
		cout << "===========" << t << "===========" <<endl;
		Utility::output_matrix(B[t],nSize,nSize);
		cout << "+++++++++++" << t+10 << "+++++++++++" <<endl;
		Utility::output_matrix(B[t+10],nSize,nSize);
		*/
		//double** B1 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices,t,nSize);
		//double** B2 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices,t+10,nSize);

		double B1 = sf->precise_digraph->rbf(x);
		double B2 = sf->precise_digraph->rbf(x+period);
		cout << B1 <<"\t" << B2 <<endl;

		/*
		if (t==1000) {
			cout << "===========" << t << "===========" <<endl;
			Utility::output_matrix(B1,nSize,nSize);
			cout << "+++++++++++" << t+10 << "+++++++++++" <<endl;
			Utility::output_matrix(B2,nSize,nSize);
		}
		*/

		/*
		for ( int i=0; i<nSize; i++) {
			for ( int j=0; j<nSize; j++) {
				if (abs(B1[i][j] +130 - B2[i][j])>0.000001) {
					equal = false;
					break;
				}
			}
			if (!equal) break;
		}
		*/

		if (abs(B1+0.13*period - B2) > 0.000001) equal = false;
		cout << "t=" << t << "=>" << equal <<endl;
		if (equal) break;
		
	}
}

/**
 *  The linear periodicity of the simplified (strongly connected) digraph with some loose edges
 */
TEST(StateflowTest, Stateflow4) {
	Digraph* digraph = DigraphExample::generateDigraph6();

	digraph->generate_strongly_connected_components();
	EXPECT_EQ(digraph->sccs.size(),1);

	cout << digraph->sccs.size() <<endl;
	for(vector<Digraph*>::iterator iter = digraph->sccs.begin(); iter != digraph->sccs.end(); iter++) {
		Digraph* subGraph = *iter;
		cout << "============Sub Graph===========" << endl;
		for (vector<Node*>::iterator nIter = subGraph->node_vec.begin(); nIter != subGraph->node_vec.end(); nIter++) {
			cout << (*nIter)->name <<endl;
		}
	}

	digraph->check_strongly_connected();
	EXPECT_EQ(digraph->strongly_connected,true);

	digraph->calculate_period_gcd();
	EXPECT_EQ(digraph->pGCD,100);

	digraph->calculate_linear_factor();
	EXPECT_EQ(digraph->linear_factor,0.13);

	digraph->tf = 10000;
	digraph->prepare_rbf_calculation(false);

	cout << "nSize=" << digraph->unit_digraph->n_size << endl;
	if(output) {
		Utility::output_matrix(digraph->unit_digraph->matrix,digraph->unit_digraph->n_size,digraph->unit_digraph->n_size);
	}

	cout << "lper=" << digraph->unit_digraph->lper << ",\tldef=" << digraph->unit_digraph->ldef <<endl;

	if (true) {
		for (int t=0; t<=1000; t+=100)
			cout << "rbf(" << t << ")=" << digraph->rbf(t) << endl;
	}


	//digraph->write_graphviz(std::cout);
}

#endif