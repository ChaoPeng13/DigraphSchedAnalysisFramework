#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"
#include "FileReader.h"
#include "FileWriter.h"
#include "StateflowExample.h"
#include "DigraphExample.h"

extern bool output;

/*
TEST(SchedulabilityAnalysisTest, DigraphSchedulable)
{
	int num = 10;
	int maxNode = 15;
	int test = 0;
	for (double total = 0.1; total<1; total+=0.1) {
		int nRBF = 0;
		int nIBF = 0;

		for (int run=0; run<10; run++) {
			Digraph** digraphs = RandomGenerator::generate_digraph_system(num, maxNode, total);
			SchedulabilityAnalysis::prepare_all_digraphs(digraphs,num);
			if (SchedulabilityAnalysis::rbf_analysis(digraphs,num,0)) nRBF++;
			//if (SchedulabilityAnalysis::ibf_analysis(digraphs,num,0)) nIBF++;
			
			// Write files
			string name = "Output\\Test"+Utility::int_to_string(test)+"Run"+Utility::int_to_string(run)+".dot";
			const char *p = name.c_str();

			FileWriter::DotFileWriter(digraphs, num, p);
			cout<<name<<"\t"<<nRBF<<"\t"<<nIBF<<endl;

			for(int i=0; i<num; i++)
				delete digraphs[i];
			delete[] digraphs;
		}
		cout<<"Utilization="<<total<<"\tnRBF="<<nRBF<<"\tnIBF="<<nIBF<<endl;
		test++;
	}
}
*/

/*
TEST(SchedulabilityAnalysisTest, DigraphSchedulable2)
{
	for (int num=0; num<=10; num++) {
		Digraph** digraphs = new Digraph*[num];
		for (int i=0; i<num; i++) {
			Digraph* digraph = DigraphExample::generateDigraph3();
			digraphs[i] = digraph;

			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();
			digraph->calculate_period_gcd();
			digraph->calculate_linear_factor();
		}

		SchedulabilityAnalysis::prepare_all_digraphs(digraphs, num);
		bool isRBF = SchedulabilityAnalysis::rbf_analysis(digraphs, num, 0);
		if (output) cout<<"NUM="<<num<<"isRBF="<<isRBF<<endl;
	
		for(int i=0; i<num; i++)
			delete digraphs[i];
		delete[] digraphs;
	}
}
*/

TEST(SchedulabilityAnalysisTest, StartTimes)
{
	Stateflow* sf0 = StateflowExample::generateStateflow0();

	sf0->calculate_gcd();
	sf0->calculate_t_gcd();
	sf0->calculate_hyperperiod();
	sf0->set_state_number();
	sf0->generate_rbf_time_instances();
	sf0->generate_rbfs();
	sf0->generate_exec_req_matrix();
	sf0->calculate_linear_factor();

	Stateflow* sf1 = StateflowExample::generateStateflow1();

	sf1->calculate_gcd();
	sf1->calculate_t_gcd();
	sf1->calculate_hyperperiod();
	sf1->set_state_number();
	sf1->generate_rbf_time_instances();
	sf1->generate_rbfs();
	sf1->generate_exec_req_matrix();
	sf1->calculate_linear_factor();

	Stateflow* sfs[2] = {sf0,sf1};
	SchedulabilityAnalysis::prepare_all_stateflows(sfs,2,0);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,2,0);
	EXPECT_EQ(isRBF, true);
}

TEST(SchedulabilityAnalysisTest, RBFArbitraryOffset)
{
	Stateflow** sfs;
	int num;
	const char* file = "Output\\Test0Run886.dot";
	FileReader::DotFileReader(sfs, num,1000,file);
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		// show rbfs
		if (output) {
			for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
				cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
				Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
			}
		}

		sf->generate_exec_req_matrix();
		// show execution request matrix
		if (output) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();
	}


	SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,0);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,0);
	bool isRBF2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,0);
	//bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,0);
	//bool isIBF2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,0);

	//cout<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
	//cout<<isRBF<<"\t"<<isRBF2<<endl; 
	EXPECT_EQ(isRBF, true);
	EXPECT_EQ(isRBF2, true);
}

TEST(SchedulabilityAnalysisTest, RBFArbitraryOffset2)
{
	// write a file
	string result = "TestResult\\result230";
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	for (int period = 230; period <=230; period++) {
		for (int wcet = 170; wcet<=170; wcet++) {
			int num = 2;
			int choice = 0;
			Stateflow** sfs = new Stateflow*[num];
			sfs[0] = StateflowExample::generateStateflow0();
			sfs[1] = StateflowExample::generateStateflow2();

			for (vector<Transition*>::iterator iter = sfs[1]->trans.begin(); iter != sfs[1]->trans.end(); iter++) {
				Transition* t = *iter;
				t->wcet = wcet;
				t->period = period;
			}
	
			for (int i=0; i<num; i++) {
				Stateflow* sf = sfs[i];
				//sf->write_graphviz(cout);
				sf->calculate_gcd();
				sf->calculate_hyperperiod();
				sf->set_state_number();
				sf->generate_rbf_time_instances();
				sf->generate_rbfs();

				sf->generate_exec_req_matrix();
				// show execution request matrix
				sf->calculate_linear_factor();
			}
			SchedulabilityAnalysis::prepare_all_stateflows(sfs, num, choice);
			bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num, choice);
			bool isRBF2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num ,choice);
			bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,choice);
			bool isIBF2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,choice);
		
			//if (wcet == 100) FileWriter::DotFileWriter(sfs[0], 1000, "DrawGraph");
			cout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			fout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			//EXPECT_EQ(isRBF, true);
			//EXPECT_EQ(isRBF2,false);
			//EXPECT_EQ(isIBF, true);
			//EXPECT_EQ(isIBF2,false);

			for (int i=0; i<num; i++)
				delete sfs[i];
			delete[] sfs;
		}
	}
	fout.close();
}

TEST(SchedulabilityAnalysisTest, RBFArbitraryOffset3)
{
	// write a file
	string result = "TestResult\\result300";
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	for (int period = 300; period <=300; period++) {
		for (int wcet = 140; wcet<=140; wcet++) {
			int num = 2;
			int choice = 0;
			Stateflow** sfs = new Stateflow*[num];
			sfs[0] = StateflowExample::generateStateflow0();
			sfs[1] = StateflowExample::generateStateflow2();

			for (vector<Transition*>::iterator iter = sfs[1]->trans.begin(); iter != sfs[1]->trans.end(); iter++) {
				Transition* t = *iter;
				t->wcet = wcet;
				t->period = period;
			}
	
			for (int i=0; i<num; i++) {
				Stateflow* sf = sfs[i];
				//sf->write_graphviz(cout);
				sf->calculate_gcd();
				sf->calculate_hyperperiod();
				sf->set_state_number();
				sf->generate_rbf_time_instances();
				sf->generate_rbfs();

				sf->generate_exec_req_matrix();
				// show execution request matrix
				sf->calculate_linear_factor();
			}
			SchedulabilityAnalysis::prepare_all_stateflows(sfs, num, choice);
			bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num, choice);
			bool isRBF2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num ,choice);
			bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,choice);
			bool isIBF2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,choice);
		
			//FileWriter::DotFileWriter(sfs[0], 1000, "DrawGraph");
			cout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			fout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			//EXPECT_EQ(isRBF, true);
			//EXPECT_EQ(isRBF2,false);
			//EXPECT_EQ(isIBF, true);
			//EXPECT_EQ(isIBF2,false);

			for (int i=0; i<num; i++)
				delete sfs[i];
			delete[] sfs;
		}
	}
	fout.close();
}

/*
TEST(SchedulabilityAnalysisTest, StateflowSchedulable)
{
	int num = 10;
	int maxState = 15;

	for (double total = 0.01; total<=0.20; total+=0.01) {
		int nRBFStatic = 0;
		int nRBFArbitrary = 0;
		int nIBFStatic = 0;
		int nIBFArbitrary = 0;

		cout<<NEG_INFINITY+200<<endl;

		for (int run=0; run<10; run++) {
			cout<<"Run-"<<run<<endl;
			Stateflow** sfs = RandomGenerator::generate_stateflow_system(num, maxState, total);
			if (SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num)) nRBFStatic++;
			if (SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num)) nRBFArbitrary++;
			//if (SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num)) nIBFStatic++;
			//if (SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num)) nIBFArbitrary++;
		}

		cout<<"Utilization="<<total<<"\tnRBFStatic="<<nRBFStatic<<"\tnRBFArbitrary="<<nRBFArbitrary<<endl;
		cout<<"            "<<total<<"\tnIBFStatic="<<nIBFStatic<<"\tnIBFArbitrary="<<nIBFArbitrary<<endl;
	}
}
*/

#endif