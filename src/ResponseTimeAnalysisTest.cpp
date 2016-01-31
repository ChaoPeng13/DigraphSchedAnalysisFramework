#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "SchedulabilityAnalysis.h"
#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"
#include "FileReader.h"
#include "FileWriter.h"
#include "StateflowExample.h"
#include "DigraphExample.h"

#include "ResponseTimeAnalysis.h"
/*
TEST(ResponseTimeAnalysisTest, ResponseTime)
{
	// write a file
	string result = "TestResult\\ResponseTime";
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

			ResponseTimeAnalysis::calculate_all_response_times(sfs,num);
		
			//FileWriter::DotFileWriter(sfs[0], 1000, "DrawGraph");
			cout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			fout<<"period="<<period<<"\t"<<"wcet="<<wcet<<"\t"<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
			//EXPECT_EQ(isRBF, true);
			//EXPECT_EQ(isRBF2,false);
			//EXPECT_EQ(isIBF, true);
			//EXPECT_EQ(isIBF2,false);

			// Write all the response times
			FileWriter::WriteResponseTimes(sfs,num,fout);
			FileWriter::WriteResponseTimes(sfs,num,cout);

			for (int i=0; i<num; i++)
				delete sfs[i];
			delete[] sfs;
		}
	}
	fout.close();
}
*/
#endif