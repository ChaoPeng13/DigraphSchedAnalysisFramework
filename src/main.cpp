/**
 *
 */
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>

#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>
#pragma comment(lib,"E:\\googletest-master\\googletest\\msvc\\gtest\\Debug\\gtestd.lib");
#endif

#ifdef VLDTEST
#include <vld.h>
#endif

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"
#include "FileWriter.h"
#include "FileReader.h"
#include "StateflowExample.h"
#include "ResponseTimeAnalysis.h"
#include "Timer.h"

void generateRandomSystemsForExactAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun) {
// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 1000;

	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=0; run<maxRun; run++) {
					try {
						Stateflow** sfs;
						while(true) {
							sfs = RandomGenerator::generate_stateflow_system_for_exact_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

							SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);
							if (SchedulabilityAnalysis::generate_request_function(sfs,num,false)) break;

							for (int i=0; i<num; i++)
								delete sfs[i];
							delete[] sfs;
						}

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
				
					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void generateRandomSystemsForApproximateAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun) {
// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;
	
	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=0; run<maxRun; run++) {
					try {
						Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
				
					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void ExactAnalysis(const char* file) {
	Stateflow** sfs;
	int num;
	//const char* file = "Input\\Stateflows\\Num3Run175.dot";
	//const char* file = "Input\\Stateflows\\Num3Run919.dot";

	// used to find an exception
	//const char* file = "Input\\Stateflows\\Num18Run534.dot";
	//const char* file = "Input\\Stateflows\\Num20Run812.dot";

	/// \brief test the original stateflow from java codes for ECRTS2012
	// used to find an exception
	// const char* file = "Origin\\OriginRun151S.dot";
	// const char* file = "Origin\\OriginRun151A.dot";
	FileReader::DotFileReader(sfs, num,1,file);
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		// show rbfs
		if (false) {
			for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
				cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
				Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
			}
		}

		sf->generate_exec_req_matrix();
		// show execution request matrix
		if (false) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();
	}


	SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,0);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,0);
	bool isRBF2 = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,0);
	bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs, num,0);
	bool isIBF2 = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,0);
	/*
	for (int i=0; i<2; i++) {
		Stateflow* sf = sfs[i];
		FileWriter::DotFileWriter(sf, 4000, "Num3Run175SF"+Utility::int_to_string(i));
	}
	*/
	cout<<isRBF<<"\t"<<isRBF2<<"\t"<<isIBF<<"\t"<<isIBF2<<endl;

	for (int i=0; i<num; i++)
		delete sfs[i];
	delete[] sfs;
}

void ExactAnalysis(string name)
{
	// write a file
	string result = "Results\\OneStateflow3";

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	string result2 = "Results\\SimpleDigraphs5";
	ofstream fout2(result2, ios::out | ios::trunc);
	if (!fout2.is_open()) {
		cerr << "Can't open "<<result2<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	string result3 = "Results\\RBF3";
	ofstream fout3(result3, ios::out | ios::trunc);
	if (!fout3.is_open()) {
		cerr << "Can't open "<<result3<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	string result4 = "Results\\sRBF3";
	ofstream fout4(result4, ios::out | ios::trunc);
	if (!fout4.is_open()) {
		cerr << "Can't open "<<result4<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int choice = 2;
	int maxState = 15;

	int nException = 0;

	Stateflow** sfs;
		
	int num;
	//string name = "Input\\TwoStateflows.dot";
	//string name = "StateflowsAndDigraphs\\Num20Util6Run929.dot";
	//string name = "StateflowsAndDigraphs2\\Num20Util7Run145.dot";
	//string name = "ExactTest\\Num20Util60Run533.dot";
	//string name = "ExactStateflowsByUtil2\\Num20Util65Scc90Run9.dot";
	const char* file = name.c_str();
	FileReader::DotFileReader(sfs, num,1,file);

	try {
		for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		// show rbfs
		if (false) {
			for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
				cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
				Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
			}
		}

		sf->generate_exec_req_matrix();
		// show execution request matrix
		if (false) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
			sf->calculate_linear_factor();
		}

		SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,choice);

		/*
		SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);

		int deadline = 0;
		for (int i=0; i<num; i++) {
			Stateflow* sfi = sfs[i];
			deadline = max(deadline,sfi->maxDeadline);
		}
		fout2 << "deadline=" << deadline <<endl;

		for (int k=1; k<=10; k++) {
			SchedulabilityAnalysis::generate_request_function(sfs[1],0,k*sfs[1]->hyperperiod);
			fout2 << k <<": "<<sfs[1]->CRF.size()<<endl;
			SchedulabilityAnalysis::output_critical_request_function(sfs[1],fout);
		}
		*/
		
		cout << "Generating critical action pairs" <<endl;
		SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);
		//SchedulabilityAnalysis::output_critical_action_pair(sfs,num,fout);
		cout << "Generating request functions" << endl;
		//SchedulabilityAnalysis::generate_request_function(sfs,num);
		//SchedulabilityAnalysis::output_critical_request_function(sfs,num,fout);
		cout << "Generating request function abstract tree" <<endl; 
		//SchedulabilityAnalysis::generate_request_function_abstract_tree(sfs,num);

		cout << "Do Exact analysis" <<endl;
		bool EXACT = true; //SchedulabilityAnalysis::exact_sched_analysis(sfs,num,choice);
		
		cout << "Do RBF analysis with static offset" <<endl;
		bool RASO = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,num,choice);
		cout << "Do IBF analysis with static offset" <<endl;
		bool IASO = false; //SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,choice);

		cout << "Do RBF analysis with arbitrary offset" <<endl;
		bool RAAO = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,choice);
		cout << "Do IBF analysis with arbitrary offset" << endl;
		bool IAAO = false; //SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,choice);

		cout << "Do RBF analysis with arbitrary offset based on simple digraph" <<endl;
		bool RAAOSD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);
		cout << "Do IBF analysis with arbitrary offset based on simple digraph" <<endl;
		bool IAAOSD = false; //SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);

		cout << "Do RBF analysis with arbitrary offset based on precise digraph" <<endl;
		bool RAAOPD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);
		cout << "Do IBF analysis with arbitrary offset based on precise digraph" <<endl;
		bool IAAOPD = false; //SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);

		cout << "Do linear RBF analysis" << endl;
		bool LURBF = SchedulabilityAnalysis::lu_rbf_sched_analysis(sfs,num,choice);
		cout << "Do linear IBF analysis" << endl;
		bool LUIBF = SchedulabilityAnalysis::lu_ibf_sched_analysis(sfs,num,choice);

		cout<<file<<endl;
		cout<<"\t"<<"=>"<<EXACT<<endl;
		cout<<"\t"<<"=>"<< RASO <<"\t"<<RAAO<<"\t"<<RAAOSD<<"\t"<<RAAOPD<<"\t"<<LURBF<<endl;
		cout<<"\t"<<"=>"<< IASO <<"\t"<<IAAO<<"\t"<<IAAOSD<<"\t"<<IAAOPD<<"\t"<<LUIBF<<endl;

		for (int i=0; i<num; i++)
			delete sfs[i];
		delete[] sfs;
	} 
	catch(bad_alloc& ba)
	{
		nException++;
		cout<<"Exception"<<endl;

		for (int i=0; i<num; i++)
			delete sfs[i];
		delete[] sfs;
	}
	fout.close();
	fout2.close();
	fout3.close();
	fout4.close();
}

void ExactAnalysis()
{
	string result = "Results\\OneStateflow2";

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	Stateflow* sf = StateflowExample::generateStateflow0();
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	sf->calculate_linear_factor();
	//sf->check_irreducible();
	sf->calculate_generialized_period();
	sf->calculate_generialized_defect();
	sf->calculate_exec_req_matrix_power(10);

	for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter++) {
		int s = *iter;
		if (s == sf->hyperperiod) continue;
		fout << "rbf[" << s <<endl;
		for ( int f=s; f<=s+sf->hyperperiod; f+=sf->gcd) {
			fout << "," << f << "]=" << sf->get_rbf(s,f) << endl;
		}
		fout << endl;
	}
	bool output = false;
	SchedulabilityAnalysis::generate_request_function(sf,0,sf->hyperperiod,output);
	SchedulabilityAnalysis::output_critical_request_function(sf,fout);
	for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter ++) {
		int stime = *iter;
		if (stime == sf->hyperperiod) continue;
		SchedulabilityAnalysis::generate_request_function_abstract_tree(sf,stime,stime+sf->hyperperiod,output);
	}
	SchedulabilityAnalysis::output_request_function_abstract_tree(sf,fout);
}

void SchedAnalysis(bool myProperty, string file, string directory,int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int startRun, int endRun) {
//#define IBF_ANALYSIS
//#define PERIODICITY_PROPERTY
	// write a file
	string result = file;
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int choice = 0;
	int maxState = 15;
	bool output = false;

	SchedulabilityAnalysis::set_zero();

	for (int num=minNum; num<=maxNum; num+=stepNum) {

		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {

			for (int scc = minScc; scc <= maxScc; scc += stepScc) {
				int nException = 0;
				int nUnSuccess = 0;

				for (int run=startRun; run<=endRun; run++) {
					Stateflow** sfs;
		
					int numStateflows;
					string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
					const char* file = name.c_str();
					FileReader::DotFileReader(sfs, numStateflows,1,file);

					try {
						Timer timer;
						timer.start();
						for (int i=0; i<numStateflows; i++) {
							Stateflow* sf = sfs[i];
							//sf->write_graphviz(cout);
							sf->calculate_gcd();
							sf->calculate_hyperperiod();
							sf->set_state_number();
							sf->generate_rbf_time_instances();
							sf->generate_rbfs();
							// show rbfs
							if (false) {
								for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
									cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
									Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
								}
							}

							// statistics the average degree for stateflows
							SchedulabilityAnalysis::avgDegree += 1.0*sf->trans.size()/sf->states.size();

							sf->generate_exec_req_matrix();
							// show execution request matrix
							if (false) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
							sf->calculate_linear_factor();
						}
						timer.end();
						SchedulabilityAnalysis::tPrepareStateflow += timer.getTime();
						
						cout << "Generating critical action pairs" <<endl;
						SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);
						if (output) SchedulabilityAnalysis::output_critical_action_pair(sfs,num,fout);
	#ifndef IBF_ANALYSIS					
						//cout << "Do RBF analysis with arbitrary offset based on simple digraph" <<endl;
						bool RAAOSD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs, num,choice);
	#endif
	#ifdef PERIODICITY_PROPERTY
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];
							sf->check_irreducible();
							if(sf->isIrred) {
								SchedulabilityAnalysis::nIrreducibleStateflows++;
								cout<<"isIrred"<<endl;
							}
						}
	#endif
						SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,2);
	#ifndef IBF_ANALYSIS					
						//cout << "Do Exact analysis" <<endl;
						bool EXACT = true;
						bool EXACTAO = true;
						if (myProperty) EXACT = SchedulabilityAnalysis::exact_sched_analysis(sfs,num,choice,false,fout);
						if (myProperty) EXACTAO = SchedulabilityAnalysis::exact_analysis_arbitrary_offset(sfs,num,choice);
						
						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						//cout << "Do RBF analysis with static offset" <<endl;
						bool RASO = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,num,choice);

						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						//cout << "Do IBF analysis with static offset" <<endl;
						bool IASO = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,choice);
						
						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						//cout << "Do RBF analysis with arbitrary offset" <<endl;
						bool RAAO = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,choice);
						
						//cout << "Do RBF analysis with arbitrary offset based on simple digraph" <<endl;
						//bool RAAOSD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);
						
						//cout << "Do RBF analysis with arbitrary offset based on precise digraph" <<endl;
						bool RAAOPD = true; // SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);
						//sfs[2]->precise_digraph->write_graphviz(cout);
						//cout << "Do linear RBF analysis" << endl;
						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						bool LURBF = SchedulabilityAnalysis::lu_rbf_sched_analysis(sfs,num,choice);
						
						//cout << "Do linear IBF analysis" << endl;
						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						bool LUIBF = SchedulabilityAnalysis::lu_ibf_sched_analysis(sfs,num,choice);

						//cout << "Do linear CSUM analysis" <<endl;
						SchedulabilityAnalysis::reset_calculating_containers(sfs,num);
						bool LUCSUM = SchedulabilityAnalysis::lu_csum_sched_analysis(sfs,num,choice);
	#endif
	#ifdef IBF_ANALYSIS
						//cout << "Do IBF analysis with arbitrary offset" << endl;
						bool IAAO = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,choice);

						//cout << "Do IBF analysis with arbitrary offset based on simple digraph" <<endl;
						bool IAAOSD = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);

						//cout << "Do IBF analysis with arbitrary offset based on precise digraph" <<endl;
						bool IAAOPD = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT) SchedulabilityAnalysis::nExactStaticOffset++;
						if (EXACTAO) SchedulabilityAnalysis::nExactArbitraryOffset++;
						
						if (RASO) SchedulabilityAnalysis::nRBFStaticOffset++;
						if (RAAO) SchedulabilityAnalysis::nRBFArbitraryOffset++;
						if (RAAOSD) SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraph++;
						if (RAAOPD) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraph++;

						if (LURBF) SchedulabilityAnalysis::nLinearUpperRBF++;
						if (LUIBF) SchedulabilityAnalysis::nLinearUpperIBF++;
						if (LUCSUM) SchedulabilityAnalysis::nLinearUpperCSUM++;

						if (IASO) SchedulabilityAnalysis::nIBFStaticOffset++;
	#endif
	#ifdef IBF_ANALYSIS
						if (IAAO) SchedulabilityAnalysis::nIBFArbitraryOffset++;
						if (IAAOSD) SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraph++;
						if (IAAOPD) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraph++;
	#endif
						//if (isRBF != isRBF2 || isRBF != isRBF3) {
	#ifndef IBF_ANALYSIS
						cout<<file<<"=>"<< EXACT << "\t" << EXACTAO <<"\t"<< RASO << "\t" << IASO <<"\t"<<RAAO<<"\t"<<RAAOSD<<"\t"<<RAAOPD
							<<"\t"<<LURBF<<"\t"<<LUIBF<<"\t"<<LUCSUM;
	#endif
	#ifdef IBF_ANALYSIS
						cout<<"\t"<<IAAO<<"\t"<<IAAOSD<<"\t"<<IAAOPD;
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT != EXACTAO) cout << "\t" << "DDD";
						if (EXACT != RASO) cout<<"\t"<<"***";
						if (RASO != RAAO) cout<<"\t"<<"###";
						if (RAAO != RAAOPD) cout<<"\t"<<"^^^";
						if (RAAOPD != RAAOSD) cout<<"\t"<<"+++";
						if (RAAOSD != LURBF) cout<<"\t"<<"@@@";
						if(LURBF != LUIBF) cout <<"\t"<<"%%%";
						if(LURBF!= LUCSUM) cout <<"\t"<<"```";
						if (IASO != RASO) cout<<"\t"<<"$$$";
	#endif
	#ifdef IBF_ANALYSIS
						
						if (IAAO != IAAOPD) cout<<"\t"<<"^x^";
						if (IAAOPD != IAAOSD) cout<<"\t"<<"+x+";
	#endif
						cout<<endl;
	#ifndef IBF_ANALYSIS
						fout<<file<<"=>"<< EXACT << "\t" <<EXACTAO <<"\t"<< RASO << "\t" << IASO <<"\t"<<RAAO<<"\t"<<RAAOSD<<"\t"<<RAAOPD
							<<"\t"<<LURBF<<"\t"<<LUIBF<<"\t"<<LUCSUM;
	#endif
	#ifdef IBF_ANALYSIS
						fout<<"\t"<<IAAO<<"\t"<<IAAOSD<<"\t"<<IAAOPD;
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT != EXACTAO) fout << "\t" << "DDD";
						if (EXACT != RASO) fout<<"\t"<<"***";
						if (RASO != RAAO) fout<<"\t"<<"###";
						if (RAAO != RAAOPD) fout<<"\t"<<"^^^";
						if (RAAOPD != RAAOSD) fout<<"\t"<<"+++";
						if (RAAOSD != LURBF) fout<<"\t"<<"@@@";
						if(LURBF != LUIBF) fout <<"\t"<<"%%%";
						if(LURBF!= LUCSUM) fout <<"\t"<<"```";
						if (IASO != RASO) fout<<"\t"<<"$$$";
	#endif
	#ifdef IBF_ANALYSIS
						if (IAAO != IAAOPD) fout<<"\t"<<"^x^";
						if (IAAOPD != IAAOSD) fout<<"\t"<<"+x+";
	#endif
						fout<<endl;

						double tUtil = 0;
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];

							SchedulabilityAnalysis::tCalDiffPeriod += sf->tDiff;
							SchedulabilityAnalysis::tCalNonPeriod += sf->tCalWithoutPeriodicity;
							SchedulabilityAnalysis::tCalPeriod += sf->tCalWithPeriodicity;

							tUtil+=sf->lfac;
						}
						cout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil<<endl;
						fout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil<<endl;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					} 
					catch(bad_alloc& ba)
					{
						nException++;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					}
				}
				fout<<"nException="<<nException<<endl;
				fout<<"nUnSuccess="<<nUnSuccess<<endl;
				SchedulabilityAnalysis::nStateflows = num;
				SchedulabilityAnalysis::totalUtilization = util;
				// reset
				SchedulabilityAnalysis::reset();
			}
		}
	}
	fout<<"Chocice="<<choice<<endl;
	
	SchedulabilityAnalysis::output_vectors(fout);

	fout.close();
}

void generateRandomSystemsForWeightedSchedAnalysis(string directory, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;
	
	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		//for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=0; run<maxRun; run++) {
					try {
						int util = rand()%(maxUtil-minUtil)+minUtil; // util \in [minUtil,maxUtil]
						Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

						string name = directory+"\\WSNum"+Utility::int_to_string(num)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
				
					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		//}
	}
}

void WeightedSchedulabilityAnalysis(bool myProperty,string directory, string file, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// #define PERIODICITY_PROPERTY
	// write a file
	string result = file;
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int choice = 2;
	int maxState = 15;
	bool output = false;

	vector<double> vec_sum_weightedSchedEXACT;
	vector<double> vec_sum_weightedSchedRASO;
	vector<double> vec_sum_weightedSchedIASO;
	vector<double> vec_sum_weightedSchedRAAO;
	vector<double> vec_sum_weightedSchedRAAOPD;
	vector<double> vec_sum_weightedSchedRAAOSD;
	vector<double> vec_sum_weightedSchedLURBF;
	vector<double> vec_sum_weightedSchedLUIBF;
	vector<double> vec_sum_weightedSchedLUCSUM;

	SchedulabilityAnalysis::set_zero();

	for (int num=minNum; num<=maxNum; num+=stepNum) {

		//for (int util = minUtil; util<=maxUtil; util+=stepUtil) {

			for (int scc = minScc; scc <= maxScc; scc += stepScc) {
				int nException = 0;

				double sum_weightedSchedEXACT = 0;
				double sum_weightedSchedRASO = 0;
				double sum_weightedSchedIASO = 0;
				double sum_weightedSchedRAAO = 0;
				double sum_weightedSchedRAAOPD = 0;
				double sum_weightedSchedRAAOSD = 0;
				double sum_weightedSchedLURBF = 0;
				double sum_weightedSchedLUIBF = 0;
				double sum_weightedSchedLUCSUM = 0;

				double totalUtil = 0;

				for (int run=0; run<maxRun; run++) {
					Stateflow** sfs;
					
					int numStateflows;
					string name = directory+"\\WSNum"+Utility::int_to_string(num)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
					const char* file = name.c_str();
					FileReader::DotFileReader(sfs, numStateflows,1,file);

					try {
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];
							//sf->write_graphviz(cout);
							sf->calculate_gcd();
							sf->calculate_hyperperiod();
							sf->set_state_number();
							sf->generate_rbf_time_instances();
							sf->generate_rbfs();
							// show rbfs
							if (false) {
								for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
									cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
									Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
								}
							}

							// statistics the average degree for stateflows
							SchedulabilityAnalysis::avgDegree += 1.0*sf->trans.size()/sf->states.size();

							sf->generate_exec_req_matrix();
							// show execution request matrix
							if (false) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
							sf->calculate_linear_factor();
						}

					
						cout << "Generating critical action pairs" <<endl;
						SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);
						if (output) SchedulabilityAnalysis::output_critical_action_pair(sfs,num,fout);
	#ifndef IBF_ANALYSIS					
						//cout << "Do RBF analysis with arbitrary offset based on simple digraph" <<endl;
						bool RAAOSD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs, num,choice);
	#endif
	#ifdef PERIODICITY_PROPERTY
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];
							sf->check_irreducible();
							if(sf->isIrred) {
								SchedulabilityAnalysis::nIrreducibleStateflows++;
								cout<<"isIrred"<<endl;
							}
						}
	#endif
						SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,choice);
	#ifndef IBF_ANALYSIS					
						//cout << "Do Exact analysis" <<endl;
						bool EXACT = true;
						if (myProperty) EXACT = SchedulabilityAnalysis::exact_sched_analysis(sfs,num,choice,false,fout);
						
						//cout << "Do RBF analysis with static offset" <<endl;
						bool RASO = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,num,choice);

						//cout << "Do IBF analysis with static offset" <<endl;
						bool IASO = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,choice);
						
						//cout << "Do RBF analysis with arbitrary offset" <<endl;
						bool RAAO = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,choice);
						
						//cout << "Do RBF analysis with arbitrary offset based on simple digraph" <<endl;
						//bool RAAOSD = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);
						
						//cout << "Do RBF analysis with arbitrary offset based on precise digraph" <<endl;
						bool RAAOPD = true; // SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);
						//sfs[2]->precise_digraph->write_graphviz(cout);
						//cout << "Do linear RBF analysis" << endl;
						bool LURBF = SchedulabilityAnalysis::lu_rbf_sched_analysis(sfs,num,choice);
						
						//cout << "Do linear IBF analysis" << endl;
						bool LUIBF = SchedulabilityAnalysis::lu_ibf_sched_analysis(sfs,num,choice);

						//cout << "Do linear CSUM analysis" <<endl;
						bool LUCSUM = SchedulabilityAnalysis::lu_csum_sched_analysis(sfs,num,choice);
	#endif
	#ifdef IBF_ANALYSIS
						//cout << "Do IBF analysis with arbitrary offset" << endl;
						bool IAAO = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,choice);

						//cout << "Do IBF analysis with arbitrary offset based on simple digraph" <<endl;
						bool IAAOSD = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,choice);

						//cout << "Do IBF analysis with arbitrary offset based on precise digraph" <<endl;
						bool IAAOPD = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,choice);
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT) SchedulabilityAnalysis::nExactStaticOffset++;
						
						if (RASO) SchedulabilityAnalysis::nRBFStaticOffset++;
						if (RAAO) SchedulabilityAnalysis::nRBFArbitraryOffset++;
						if (RAAOSD) SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraph++;
						if (RAAOPD) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraph++;

						if (LURBF) SchedulabilityAnalysis::nLinearUpperRBF++;
						if (LUIBF) SchedulabilityAnalysis::nLinearUpperIBF++;
						if (LUCSUM) SchedulabilityAnalysis::nLinearUpperCSUM++;

						if (IASO) SchedulabilityAnalysis::nIBFStaticOffset++;
	#endif
	#ifdef IBF_ANALYSIS
						if (IAAO) SchedulabilityAnalysis::nIBFArbitraryOffset++;
						if (IAAOSD) SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraph++;
						if (IAAOPD) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraph++;
	#endif
						//if (isRBF != isRBF2 || isRBF != isRBF3) {
	#ifndef IBF_ANALYSIS
						cout<<file<<"=>"<< EXACT<<"\t"<< RASO << "\t" << IASO <<"\t"<<RAAO<<"\t"<<RAAOSD<<"\t"<<RAAOPD
							<<"\t"<<LURBF<<"\t"<<LUIBF<<"\t"<<LUCSUM;
	#endif
	#ifdef IBF_ANALYSIS
						cout<<"\t"<<IAAO<<"\t"<<IAAOSD<<"\t"<<IAAOPD;
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT != RASO) cout<<"\t"<<"***";
						if (RASO != RAAO) cout<<"\t"<<"###";
						if (RAAO != RAAOPD) cout<<"\t"<<"^^^";
						if (RAAOPD != RAAOSD) cout<<"\t"<<"+++";
						if (RAAOSD != LURBF) cout<<"\t"<<"@@@";
						if(LURBF != LUIBF) cout <<"\t"<<"%%%";
						if(LURBF!= LUCSUM) cout <<"\t"<<"```";
						if (IASO != RASO) cout<<"\t"<<"$$$";
	#endif
	#ifdef IBF_ANALYSIS
						
						if (IAAO != IAAOPD) cout<<"\t"<<"^x^";
						if (IAAOPD != IAAOSD) cout<<"\t"<<"+x+";
	#endif
						cout<<endl;
	#ifndef IBF_ANALYSIS
						fout<<file<<"=>"<< EXACT<<"\t"<< RASO << "\t" << IASO <<"\t"<<RAAO<<"\t"<<RAAOSD<<"\t"<<RAAOPD
							<<"\t"<<LURBF<<"\t"<<LUIBF<<"\t"<<LUCSUM;
	#endif
	#ifdef IBF_ANALYSIS
						fout<<"\t"<<IAAO<<"\t"<<IAAOSD<<"\t"<<IAAOPD;
	#endif
	#ifndef IBF_ANALYSIS
						if (EXACT != RASO) fout<<"\t"<<"***";
						if (RASO != RAAO) fout<<"\t"<<"###";
						if (RAAO != RAAOPD) fout<<"\t"<<"^^^";
						if (RAAOPD != RAAOSD) fout<<"\t"<<"+++";
						if (RAAOSD != LURBF) fout<<"\t"<<"@@@";
						if(LURBF != LUIBF) fout <<"\t"<<"%%%";
						if(LURBF!= LUCSUM) fout <<"\t"<<"```";
						if (IASO != RASO) fout<<"\t"<<"$$$";
	#endif
	#ifdef IBF_ANALYSIS
						if (IAAO != IAAOPD) fout<<"\t"<<"^x^";
						if (IAAOPD != IAAOSD) fout<<"\t"<<"+x+";
	#endif
						fout<<endl;

						double tUtil = 0;
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];

							SchedulabilityAnalysis::tCalDiffPeriod += sf->tDiff;
							SchedulabilityAnalysis::tCalNonPeriod += sf->tCalWithoutPeriodicity;
							SchedulabilityAnalysis::tCalPeriod += sf->tCalWithPeriodicity;

							tUtil+=sf->lfac;
						}

						totalUtil += tUtil;

						sum_weightedSchedEXACT += tUtil * EXACT;
						sum_weightedSchedRASO += tUtil * RASO;
						sum_weightedSchedIASO += tUtil * IASO;
						sum_weightedSchedRAAO += tUtil * RAAO;
						sum_weightedSchedRAAOPD += tUtil * RAAOPD;
						sum_weightedSchedRAAOSD += tUtil * RAAOSD;
						sum_weightedSchedLURBF += tUtil * LURBF;
						sum_weightedSchedLUIBF += tUtil * LUIBF;
						sum_weightedSchedLUCSUM += tUtil * LUCSUM;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					} 
					catch(bad_alloc& ba)
					{
						nException++;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					}
				}
				fout<<"nException="<<nException<<endl;
				SchedulabilityAnalysis::nStateflows = num;

				vec_sum_weightedSchedEXACT.push_back(sum_weightedSchedEXACT/totalUtil);
				vec_sum_weightedSchedRASO.push_back(sum_weightedSchedRASO/totalUtil);
				vec_sum_weightedSchedIASO.push_back(sum_weightedSchedIASO/totalUtil);
				vec_sum_weightedSchedRAAO.push_back(sum_weightedSchedRAAO/totalUtil);
				vec_sum_weightedSchedRAAOPD.push_back(sum_weightedSchedRAAOPD/totalUtil);
				vec_sum_weightedSchedRAAOSD.push_back(sum_weightedSchedRAAOSD/totalUtil);
				vec_sum_weightedSchedLURBF.push_back(sum_weightedSchedLURBF/totalUtil);
				vec_sum_weightedSchedLUIBF.push_back(sum_weightedSchedLUIBF/totalUtil);
				vec_sum_weightedSchedLUCSUM.push_back(sum_weightedSchedLUCSUM/totalUtil);

				//SchedulabilityAnalysis::totalUtilization = util;
				// reset
				SchedulabilityAnalysis::reset();
			}
		//}
	}
	fout<<"Chocice="<<choice<<endl;
	
	SchedulabilityAnalysis::output_vectors(fout);

	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedEXACT",vec_sum_weightedSchedEXACT);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedRASO",vec_sum_weightedSchedRASO);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedIASO",vec_sum_weightedSchedIASO);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedRAAO",vec_sum_weightedSchedRAAO);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedRAAOPD",vec_sum_weightedSchedRAAOPD);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedRAAOSD",vec_sum_weightedSchedRAAOSD);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedLURBF",vec_sum_weightedSchedLURBF);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedLUIBF",vec_sum_weightedSchedLUIBF);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSchedLUCSUM",vec_sum_weightedSchedLUCSUM);

	fout.close();
}

int main(int argc, char* argv[])
{
#ifdef GOOGLETEST
	// Run Google test
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
#endif
	string index = "20160127B2";
	string result = "Results\\";
	int period_choice = RandomGenerator::ExactAnalysis3;
	//int period_choice = RandomGenerator::ApproximateAnalysis1;
	// myProperty = true, do exact analysis
	// otherwise, do approximate analysis
	bool myProperty = true; 
	int times = 1000;

	if (true) {
		if (myProperty) {
			generateRandomSystemsForExactAnalysis("ExactAnalysisUtil"+index,period_choice,10,10,1,5,95,5,90,90,10,times);
			SchedAnalysis(myProperty,result+"ExactAnalysisUtil"+index,"ExactAnalysisUtil"+index,10,10,1,5,95,5,90,90,10,0,times-1);
		}
		else {
			//generateRandomSystemsForApproximateAnalysis("ApproximateAnalysisUtil"+index,period_choice,20,20,1,5,95,5,90,90,10,times);
			SchedAnalysis(myProperty,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,95,5,90,90,10,0,times-1);
			//SchedAnalysis(myProperty,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,35,35,5,90,90,10,5,5);
		}
	} else {
		if (myProperty) {
			//generateRandomSystemsForExactAnalysis("ExactAnalysisNum"+index,period_choice,2,20,1,50,50,5,90,90,10,times);
			SchedAnalysis(myProperty,result+"ExactAnalysisNum"+index,"ExactAnalysisNum"+index,2,20,1,50,50,5,90,90,10,0,times-1);
		}
		else {
			//generateRandomSystemsForApproximateAnalysis("ApproximateAnalysisNum"+index,period_choice,5,60,5,50,50,5,90,90,10,times);
			SchedAnalysis(myProperty,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,50,50,5,90,90,10,0,times-1);
		}
	}
	
	/*
	generateRandomSystemsForWeightedSchedAnalysis("WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
	WeightedSchedulabilityAnalysis(myProperty,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
	*/
	return 0;
}