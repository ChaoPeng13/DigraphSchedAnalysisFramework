/**
 *
 */

//#include <gtest/gtest.h>
#include <vld.h>

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"
#include "FileWriter.h"
#include "FileReader.h"

#define RUN_TEST false

void test(){
	Stateflow** sfs;
	int num;
	//const char* file = "Output\\Test0Run886.dot";
	const char* file = "Output\\Test0Run69.dot";
	FileReader::DotFileReader(sfs, num,1,file);
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

	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num, 2);
	bool isRBF2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num ,2);
	bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,2);
	bool isIBF2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,2);

	cout<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;
	for (int i=0; i<num; i++)
		delete sfs[i];
	delete[] sfs;
}

int main(int argc, char* argv[])
{
	//if (RUN_TEST) {
		// Run Google test
	//	testing::InitGoogleTest(&argc, argv);
	//	return RUN_ALL_TESTS();
	//}

	//test();
#if 1
	int choice = 0;
	//int num = 20;
	int maxState = 15;
	int test =0;
	for (int num=20; num<=20; num++) {
		for (double total = 0.07; total<0.08; total+=0.01) {
			int nRBFStatic = 0;
			int nRBFArbitrary = 0;
			int nIBFStatic = 0;
			int nIBFArbitrary = 0;

			for (int run=0; run<1000; run++) {
				//cout<<"Run-"<<run<<endl;
				try {
					Stateflow** sfs = RandomGenerator::generate_stateflow_system(num, maxState, total);
					//if (run == 886) {
					//cout<<"Now rbf..."<<endl;
					/* test all busy period length
					bool isRBFc0 = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,0);
					bool isRBFc1 = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,1);
					bool isRBFc2 = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,2);
					if (isRBFc0 != isRBFc1 || isRBFc1 != isRBFc2) {
						cout<<"Run="<<run<<"\t"<<"rbf analysis with static offsets"<<endl;
						cout<<isRBFc0<<"\t"<<isRBFc1<<"\t"<<isRBFc2<<endl;
					}
				
					bool isRBF2c0 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num,0);
					bool isRBF2c1 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num,1);
					bool isRBF2c2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num,2);
					if (isRBF2c0 != isRBF2c1 || isRBF2c1 != isRBF2c2) {
						cout<<"Run="<<run<<"\t"<<"rbf analysis with arbitary offsets"<<endl;
						cout<<isRBF2c0<<"\t"<<isRBF2c1<<"\t"<<isRBF2c2<<endl;
					}

					bool isIBFc0 = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,0);
					bool isIBFc1 = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,1);
					bool isIBFc2 = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,2);
					if (isIBFc0 != isIBFc1 || isIBFc1 != isIBFc2) {
						cout<<"Run="<<run<<"\t"<<"ibf analysis with static offsets"<<endl;
						cout<<isIBFc0<<"\t"<<isIBFc1<<"\t"<<isIBFc2<<endl;
					}
				
					bool isIBF2c0 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs, num,0);
					bool isIBF2c1 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs, num,1);
					bool isIBF2c2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs, num,2);
					if (isIBF2c0 != isIBF2c1 || isIBF2c1 != isIBF2c2) {
						cout<<"Run="<<run<<"\t"<<"ibf analysis with arbitary offsets"<<endl;
						cout<<isIBF2c0<<"\t"<<isIBF2c1<<"\t"<<isIBF2c2<<endl;
					}
					*/

					/*
					if (SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num,2)) nRBFArbitrary++;
					//cout<<"Now ibf..."<<endl;
					bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,2);
					if (isIBF) nIBFStatic++;
					if (SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,2)) nIBFArbitrary++;
				
					if (isRBF != isIBF) {
						string name = "Output\\Test"+Utility::int_to_string(test)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						FileWriter::DotFileWriter(sfs, num, p);
						cout<<name<<"\t"<<isRBF<<"\t"<<isIBF<<endl;
					}
					*/

					if (SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,choice)) nRBFStatic++;
					if (SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs, num, choice)) nRBFArbitrary++;
					if (SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,choice)) nIBFStatic++;
					if (SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num, choice)) nIBFArbitrary++;

					for (int i=0; i<num; i++)
						delete sfs[i];
					delete[] sfs;
				
				}
				catch(bad_alloc) // catch bad_alloc exception
				{
					cout<<"Run="<<run<<"\tException"<<endl;
					cout<<"Utilization="<<total<<"\tnRBFStatic="<<nRBFStatic<<"\tnRBFArbitrary="<<nRBFArbitrary<<endl;
					cout<<"            "<<total<<"\tnIBFStatic="<<nIBFStatic<<"\tnIBFArbitrary="<<nIBFArbitrary<<endl;
					return 0;
				}
			}
			test++;
			cout<<"NUM="<<num<<endl;
			cout<<"Utilization="<<total<<"\tnRBFStatic="<<nRBFStatic<<"\tnRBFArbitrary="<<nRBFArbitrary<<endl;
			cout<<"            "<<total<<"\tnIBFStatic="<<nIBFStatic<<"\tnIBFArbitrary="<<nIBFArbitrary<<endl;

			cout<<"Choice="<<choice<<endl;
			cout<<"bpRatio0="<<SchedulabilityAnalysis::bpRatio0/1000<<"\tbpRatio1="<<SchedulabilityAnalysis::bpRatio1/1000<<endl;
			cout<<"tCalCSum="<<SchedulabilityAnalysis::tCalCSum<<"\t"
				<<"tCalLinearBounds="<<SchedulabilityAnalysis::tCalLinearBounds<<endl;
			cout<<"tCalTF0="<<SchedulabilityAnalysis::tCalTF0<<"\t"
				<<"tCalTF1="<<SchedulabilityAnalysis::tCalTF1<<"\t"
				<<"tCalTF2="<<SchedulabilityAnalysis::tCalTF2<<endl;
			cout<<"tRBFStaticOffset="<<SchedulabilityAnalysis::tRBFStaticOffset<<"\t"
				<<"tRBFArbitraryOffset="<<SchedulabilityAnalysis::tRBFArbitraryOffset<<"\t"
				<<"tIBFStaticOffset="<<SchedulabilityAnalysis::tIBFStaticOffset<<"\t"
				<<"tIBFArbitraryOffset="<<SchedulabilityAnalysis::tIBFArbitraryOffset<<endl;

			// reset
			SchedulabilityAnalysis::bpRatio0 = 0;
			SchedulabilityAnalysis::bpRatio1 = 0;

			SchedulabilityAnalysis::tCalCSum = 0;
			SchedulabilityAnalysis::tCalLinearBounds = 0;
			SchedulabilityAnalysis::tCalTF0 = 0;
			SchedulabilityAnalysis::tCalTF1 = 0;
			SchedulabilityAnalysis::tCalTF2 = 0;

			SchedulabilityAnalysis::tRBFStaticOffset = 0;
			SchedulabilityAnalysis::tRBFArbitraryOffset = 0;
			SchedulabilityAnalysis::tIBFStaticOffset = 0;
			SchedulabilityAnalysis::tIBFArbitraryOffset = 0;
		}
	}
#endif

	return 0;
}