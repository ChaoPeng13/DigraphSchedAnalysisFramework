/**
 *
 */

#include <gtest/gtest.h>

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"

#define RUN_TEST false

int main(int argc, char* argv[])
{
	if (RUN_TEST) {
		// Run Google test
		testing::InitGoogleTest(&argc, argv);
		return RUN_ALL_TESTS();
	}

	int num = 10;
	int maxState = 15;

	for (double total = 0.1; total<=0.20; total+=0.01) {
		int nRBFStatic = 0;
		int nRBFArbitrary = 0;
		int nIBFStatic = 0;
		int nIBFArbitrary = 0;

		for (int run=0; run<100; run++) {
			//cout<<"Run-"<<run<<endl;
			Stateflow** sfs = RandomGenerator::generate_stateflow_system(num, maxState, total);
			if (SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num)) nRBFStatic++;
			//if (SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num)) nRBFArbitrary++;
			//if (SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num)) nIBFStatic++;
			//if (SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num)) nIBFArbitrary++;
			//delete[] sfs;
		}

		cout<<"Utilization="<<total<<"\tnRBFStatic="<<nRBFStatic<<"\tnRBFArbitrary="<<nRBFArbitrary<<endl;
		cout<<"            "<<total<<"\tnIBFStatic="<<nIBFStatic<<"\tnIBFArbitrary="<<nIBFArbitrary<<endl;
	}

	return 0;
}