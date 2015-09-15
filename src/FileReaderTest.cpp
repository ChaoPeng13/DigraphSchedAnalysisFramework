#if 0
#include <gtest/gtest.h>

#include "FileReader.h"
#include "SchedulabilityAnalysis.h"

extern bool output;

TEST(FileReaderTest, Test0)
{
	Stateflow** sfs;
	int num;
	const char* file = "Output\\Run0.dot";
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

	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num);
	bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num);

	cout<<isRBF<<"\t"<<isIBF<<endl;
}
#endif