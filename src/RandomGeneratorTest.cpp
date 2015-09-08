#include <gtest/gtest.h>

#include "RandomGenerator.h"

extern bool output;

TEST(RandomGeneratorTest, RandomDigraph)
{
	//Digraph* digraph = RandomGenerator::generate_one_digraph(0,RandomGenerator::DIGRAPH_SCALE,5,5);

	//digraph->write_graphviz(cout);
	for (int run=0; run<1; run++) {

	Digraph** digraphs = RandomGenerator::generate_digraph_system(10,15,1.0);

	for (int i=0; i<10; i++) {
		Digraph* digraph = digraphs[i];

		//digraph->write_graphviz(cout);
		cout<<"Digraph-"<<i<<": utilization=" << digraph->linear_factor <<endl;
	}
	}
}


TEST(RandomGeneratorTest, RandomStateflow)
{
	EXPECT_EQ(50000,Utility::math_lcm(50000,50000));
	int count = 0;

	for (int run=0; run<0; run++) {
		Stateflow** sfs = RandomGenerator::generate_stateflow_system(10,15,0.2);

		for (int i=0; i<10; i++) {
			Stateflow* sf = sfs[i];
			sf->calculate_t_gcd();
		//cout<<sf->gcd<<"\t";
			sf->generate_simple_digraph();
		//cout<<sf->lfac<<"\t";
			sf->generate_precise_digraph();
		//cout<<sf->precise_digraph->linear_factor<<"\t";
		
			EXPECT_EQ(int(sf->lfac*10000),int(sf->precise_digraph->linear_factor*10000));
		//sf->precise_digraph->write_graphviz(cout);
		//sf->write_graphviz(cout);
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		//sf->exec_digraph->write_graphviz(cout);
			//sf->calculate_linear_factor();
		//if (i==1) break;
			sf->check_irreducible();

			if (sf->isIrred) count++;

			sf->calculate_generialized_period();

			sf->tf0 = 100;
			sf->calculate_generialized_defect();
			cout<<"Index "<<i<<": "<<"Irreducible="<<sf->isIrred<<"\tlfac="<<sf->lfac<<"\tgper="<<sf->gper<<"\tgdef="<<sf->gdef<<endl;
			if (sf->gdef==-1) 
				Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);

			sf->calculate_exec_req_matrix_power(10);

			if (output) {
				for (int i=0; i<2000; i+=100) {
					for (int j=0; j<2000; j+=100)
						cout<<"rbf["<<+i<<","<<j<<")="<<sf->get_rbf(i,j)<<endl;
				}
			}

			for (int s=0; s<10000; s+=sf->gcd) 	for (int f=0; f<10000; f+=sf->gcd) {
				if (sf->get_rbf(s,f)!=sf->get_ibf(s,f))
					cout<<"Error goes here:"<<s<<","<<f<<endl;
				EXPECT_EQ(sf->get_rbf(s,f),sf->get_ibf(s,f));	
				EXPECT_EQ(sf->get_rbf(s,f),sf->get_dbf(s,f));
			}

		// test dbf(t)
			for (int t=0; t<10000; t+=sf->gcd)
				EXPECT_EQ(sf->get_rbf(t),sf->get_dbf(t));

			for (int i=0; i<1000; i+=10) for (int j=0; j<1000; j+=10) 
				EXPECT_LE(sf->get_ibf(i,j),sf->get_rbf(i,j));

			// test ibf(t)
			for (int t=0; t<10000; t+=sf->gcd)
				EXPECT_EQ(sf->get_rbf(t),sf->get_ibf(t));

			for (int t=0; t<1000; t+=10)
				EXPECT_GE(sf->get_rbf(t),sf->get_ibf(t));
		}

		//sf->write_graphviz(cout);
	}
	cout<<"Irreducible number = "<<count<<endl;
}