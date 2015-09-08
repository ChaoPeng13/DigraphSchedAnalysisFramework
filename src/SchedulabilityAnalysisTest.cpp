#include <gtest/gtest.h>

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"

TEST(SchedulabilityAnalysisTest, DigraphSchedulable)
{
	int num = 10;
	int maxNode = 15;

	for (double total = 0.1; total<0; total+=0.1) {
		int nRBF = 0;
		int nIBF = 0;

		for (int run=0; run<1; run++) {
			Digraph** digraphs = RandomGenerator::generate_digraph_system(num, maxNode, total);
			if (SchedulabilityAnalysis::rbf_analysis(digraphs,num)) nRBF++;
			if (SchedulabilityAnalysis::ibf_analysis(digraphs,num)) nIBF++;
		}
		cout<<"Utilization="<<total<<"\tnRBF="<<nRBF<<"\tnIBF"<<nIBF<<endl;
	}
}

TEST(SchedulabilityAnalysisTest, StateflowSchedulable)
{
	int num = 10;
	int maxState = 15;

	for (double total = 0.01; total<=0.20; total+=0.01) {
		int nRBFStatic = 0;
		int nRBFArbitrary = 0;
		int nIBFStatic = 0;
		int nIBFArbitrary = 0;

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