#include "SchedulabilityAnalysis.h"
#include "Timer.h"

#include <algorithm>

double SchedulabilityAnalysis::tfRatio0 = 0;
double SchedulabilityAnalysis::tfRatio1 = 0;

int SchedulabilityAnalysis::nIrreducibleDigraphs =0;

double SchedulabilityAnalysis::tDigraphCalSum = 0;
double SchedulabilityAnalysis::tDigraphCalLinearBounds = 0;
double SchedulabilityAnalysis::tDigraphCalTF0 = 0;
double SchedulabilityAnalysis::tDigraphCalTF1 = 0;
double SchedulabilityAnalysis::tDigraphCalTF2 = 0;
double SchedulabilityAnalysis::tDigraphRBF = 0;
double SchedulabilityAnalysis::tDigraphIBF = 0;

int SchedulabilityAnalysis::nStateflows = 0;
double SchedulabilityAnalysis::totalUtilization = 0;
double SchedulabilityAnalysis::avgDegree = 0;

int SchedulabilityAnalysis::nExactStaticOffset = 0;
int SchedulabilityAnalysis::nExactArbitraryOffset = 0;
int SchedulabilityAnalysis::nRBFStaticOffset = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffset = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraph = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraph = 0;
int SchedulabilityAnalysis::nIBFStaticOffset = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffset = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraph = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraph = 0;
int SchedulabilityAnalysis::nLinearUpperRBF = 0;
int SchedulabilityAnalysis::nLinearUpperIBF = 0;
int SchedulabilityAnalysis::nLinearUpperCSUM = 0;

double SchedulabilityAnalysis::bpRatio0 = 0;
double SchedulabilityAnalysis::bpRatio1 = 0;

int SchedulabilityAnalysis::nIrreducibleStateflows = 0;

double SchedulabilityAnalysis::tPrepareStateflow = 0;

double SchedulabilityAnalysis::tCalCSum = 0;
double SchedulabilityAnalysis::tCalLinearBounds = 0;
double SchedulabilityAnalysis::tCalTF0 = 0;
double SchedulabilityAnalysis::tCalTF1 = 0;
double SchedulabilityAnalysis::tCalTF2 = 0;

double SchedulabilityAnalysis::tCalNonPeriod = 0;
double SchedulabilityAnalysis::tCalPeriod = 0;
double SchedulabilityAnalysis::tCalDiffPeriod = 0;

double SchedulabilityAnalysis::tCalRequestExecutionMatrixWithoutPeriodicityProperty = 0;
double SchedulabilityAnalysis::tCalRequestExecutionMatrixWithPeriodicityProperty = 0;

double SchedulabilityAnalysis::tGenerateCriticalActionPairs = 0;
double SchedulabilityAnalysis::tGenerateRequestFunctions = 0;
double SchedulabilityAnalysis::tGenerateRequestFunctionAbstractTree = 0;

double SchedulabilityAnalysis::tExactStaticOffset = 0;
double SchedulabilityAnalysis::tExactArbitraryOffset = 0;
double SchedulabilityAnalysis::tRBFStaticOffset = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffset = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffsetBySimpleDigraph = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffsetByPreciseDigraph = 0;
double SchedulabilityAnalysis::tIBFStaticOffset = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffset = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetBySimpleDigraph = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetByPreciseDigraph = 0;
double SchedulabilityAnalysis::tLinearUpperRBF = 0;
double SchedulabilityAnalysis::tLinearUpperIBF = 0;
double SchedulabilityAnalysis::tLinearUpperCSUM = 0;

vector<double> SchedulabilityAnalysis::vec_tfRatio0;
vector<double> SchedulabilityAnalysis::vec_tfRatio1;

vector<int> SchedulabilityAnalysis::vec_nIrreducibleDigraphs;

vector<double> SchedulabilityAnalysis::vec_tDigraphCalSum;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalLinearBounds;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF0;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF1;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF2;
vector<double> SchedulabilityAnalysis::vec_tDigraphRBF;
vector<double> SchedulabilityAnalysis::vec_tDigraphIBF;

vector<int> SchedulabilityAnalysis::vec_nStateflows;
vector<double> SchedulabilityAnalysis::vec_totalUtilization;
vector<double> SchedulabilityAnalysis::vec_avgDegree;

vector<int> SchedulabilityAnalysis::vec_nExactStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nExactArbitraryOffset;
vector<int> SchedulabilityAnalysis::vec_nRBFStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffset;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetBySimpleDigraph;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetByPreciseDigraph;
vector<int> SchedulabilityAnalysis::vec_nIBFStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffset;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetBySimpleDigraph;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetByPreciseDigraph;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperRBF;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperIBF;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperCSUM;

vector<double> SchedulabilityAnalysis::vec_bpRatio0;
vector<double> SchedulabilityAnalysis::vec_bpRatio1;

vector<int> SchedulabilityAnalysis::vec_nIrreducibleStateflows;

vector<double> SchedulabilityAnalysis::vec_tPrepareStateflow;

vector<double> SchedulabilityAnalysis::vec_tCalCSum;
vector<double> SchedulabilityAnalysis::vec_tCalLinearBounds;
vector<double> SchedulabilityAnalysis::vec_tCalTF0;
vector<double> SchedulabilityAnalysis::vec_tCalTF1;
vector<double> SchedulabilityAnalysis::vec_tCalTF2;

vector<double> SchedulabilityAnalysis::vec_tCalNonPeriod;
vector<double> SchedulabilityAnalysis::vec_tCalPeriod;
vector<double> SchedulabilityAnalysis::vec_tCalDiffPeriod;

vector<double> SchedulabilityAnalysis::vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty;
vector<double> SchedulabilityAnalysis::vec_tCalRequestExecutionMatrixWithPeriodicityProperty;

vector<double> SchedulabilityAnalysis::vec_tGenerateCriticalActionPairs;
vector<double> SchedulabilityAnalysis::vec_tGenerateRequestFunctions;
vector<double> SchedulabilityAnalysis::vec_tGenerateRequestFunctionAbstractTree;

vector<double> SchedulabilityAnalysis::vec_tExactStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tExactArbitraryOffset;
vector<double> SchedulabilityAnalysis::vec_tRBFStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffset;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetBySimpleDigraph;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetByPreciseDigraph;
vector<double> SchedulabilityAnalysis::vec_tIBFStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffset;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetBySimpleDigraph;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetByPreciseDigraph;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperRBF;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperIBF;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperCSUM;

void SchedulabilityAnalysis::save_results() {
	vec_tfRatio0.push_back(tfRatio0);
	vec_tfRatio1.push_back(tfRatio1);

	vec_nIrreducibleDigraphs.push_back(nIrreducibleDigraphs);

	vec_tDigraphCalSum.push_back(tDigraphCalSum);
	vec_tDigraphCalLinearBounds.push_back(tDigraphCalLinearBounds);
	vec_tDigraphCalTF0.push_back(tDigraphCalTF0);
	vec_tDigraphCalTF1.push_back(tDigraphCalTF1);
	vec_tDigraphCalTF2.push_back(tDigraphCalTF2);
	vec_tDigraphRBF.push_back(tDigraphRBF);
	vec_tDigraphIBF.push_back(tDigraphIBF);

	vec_nStateflows.push_back(nStateflows);
	vec_totalUtilization.push_back(totalUtilization);
	vec_avgDegree.push_back(avgDegree);
	
	vec_nExactStaticOffset.push_back(nExactStaticOffset);
	vec_nExactArbitraryOffset.push_back(nExactArbitraryOffset);
	vec_nRBFStaticOffset.push_back(nRBFStaticOffset);
	vec_nRBFArbitraryOffset.push_back(nRBFArbitraryOffset);
	vec_nRBFArbitraryOffsetBySimpleDigraph.push_back(nRBFArbitraryOffsetBySimpleDigraph);
	vec_nRBFArbitraryOffsetByPreciseDigraph.push_back(nRBFArbitraryOffsetByPreciseDigraph);
	vec_nIBFStaticOffset.push_back(nIBFStaticOffset);
	vec_nIBFArbitraryOffset.push_back(nIBFArbitraryOffset);
	vec_nIBFArbitraryOffsetBySimpleDigraph.push_back(nIBFArbitraryOffsetBySimpleDigraph);
	vec_nIBFArbitraryOffsetByPreciseDigraph.push_back(nIBFArbitraryOffsetByPreciseDigraph);
	vec_nLinearUpperRBF.push_back(nLinearUpperRBF);
	vec_nLinearUpperIBF.push_back(nLinearUpperIBF);
	vec_nLinearUpperCSUM.push_back(nLinearUpperCSUM);

	vec_bpRatio0.push_back(bpRatio0);
	vec_bpRatio1.push_back(bpRatio1);

	vec_nIrreducibleStateflows.push_back(nIrreducibleStateflows);

	vec_tPrepareStateflow.push_back(tPrepareStateflow);

	vec_tCalCSum.push_back(tCalCSum);
	vec_tCalLinearBounds.push_back(tCalLinearBounds);
	vec_tCalTF0.push_back(tCalTF0);
	vec_tCalTF1.push_back(tCalTF1);
	vec_tCalTF2.push_back(tCalTF2);

	vec_tCalNonPeriod.push_back(tCalNonPeriod);
	vec_tCalPeriod.push_back(tCalPeriod);
	vec_tCalDiffPeriod.push_back(tCalDiffPeriod);

	vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty.push_back(tCalRequestExecutionMatrixWithoutPeriodicityProperty);
	vec_tCalRequestExecutionMatrixWithPeriodicityProperty.push_back(tCalRequestExecutionMatrixWithPeriodicityProperty);

	vec_tGenerateCriticalActionPairs.push_back(tGenerateCriticalActionPairs);
	vec_tGenerateRequestFunctions.push_back(tGenerateRequestFunctions);
	vec_tGenerateRequestFunctionAbstractTree.push_back(tGenerateRequestFunctionAbstractTree);
	
	vec_tExactStaticOffset.push_back(tExactStaticOffset);
	vec_tExactArbitraryOffset.push_back(tExactArbitraryOffset);
	vec_tRBFStaticOffset.push_back(tRBFStaticOffset);
	vec_tRBFArbitraryOffset.push_back(tRBFArbitraryOffset);
	vec_tRBFArbitraryOffsetBySimpleDigraph.push_back(tRBFArbitraryOffsetBySimpleDigraph);
	vec_tRBFArbitraryOffsetByPreciseDigraph.push_back(tRBFArbitraryOffsetByPreciseDigraph);
	vec_tIBFStaticOffset.push_back(tIBFStaticOffset);
	vec_tIBFArbitraryOffset.push_back(tIBFArbitraryOffset);
	vec_tIBFArbitraryOffsetBySimpleDigraph.push_back(tIBFArbitraryOffsetBySimpleDigraph);
	vec_tIBFArbitraryOffsetByPreciseDigraph.push_back(tIBFArbitraryOffsetByPreciseDigraph);
	vec_tLinearUpperRBF.push_back(tLinearUpperRBF);
	vec_tLinearUpperIBF.push_back(tLinearUpperIBF);
	vec_tLinearUpperCSUM.push_back(tLinearUpperCSUM);
}

void SchedulabilityAnalysis::set_zero() {
	// reset all the parameters (=0)
	tfRatio0 = 0;
	tfRatio1 = 0;

	nIrreducibleDigraphs = 0;

	tDigraphCalSum = 0;
	tDigraphCalLinearBounds = 0;
	tDigraphCalTF0 = 0;
	tDigraphCalTF1 = 0;
	tDigraphCalTF2 = 0;
	tDigraphRBF = 0;
	tDigraphIBF = 0;

	nStateflows = 0;
	totalUtilization = 0;
	avgDegree = 0;

	nExactStaticOffset = 0;
	nExactArbitraryOffset = 0;
	nRBFStaticOffset = 0;
	nRBFArbitraryOffset = 0;
	nRBFArbitraryOffsetBySimpleDigraph = 0;
	nRBFArbitraryOffsetByPreciseDigraph = 0;
	nIBFStaticOffset = 0;
	nIBFArbitraryOffset = 0;
	nIBFArbitraryOffsetBySimpleDigraph = 0;
	nIBFArbitraryOffsetByPreciseDigraph = 0;
	nLinearUpperRBF = 0;
	nLinearUpperIBF = 0;
	nLinearUpperCSUM = 0;

	bpRatio0 = 0;
	bpRatio1 = 0;

	nIrreducibleStateflows = 0;

	tPrepareStateflow = 0;

	tCalCSum = 0;
	tCalLinearBounds = 0;
	tCalTF0 = 0;
	tCalTF1 = 0;
	tCalTF2 = 0;

	tCalNonPeriod = 0;
	tCalPeriod = 0;
	tCalDiffPeriod = 0;

	tCalRequestExecutionMatrixWithoutPeriodicityProperty = 0;
	tCalRequestExecutionMatrixWithPeriodicityProperty = 0;

	tGenerateCriticalActionPairs = 0;
	tGenerateRequestFunctions = 0;
	tGenerateRequestFunctionAbstractTree = 0;

	tExactStaticOffset = 0;
	tExactArbitraryOffset = 0;
	tRBFStaticOffset = 0;
	tRBFArbitraryOffset = 0;
	tRBFArbitraryOffsetBySimpleDigraph = 0;
	tRBFArbitraryOffsetByPreciseDigraph = 0;
	tIBFStaticOffset = 0;
	tIBFArbitraryOffset = 0;
	tIBFArbitraryOffsetBySimpleDigraph = 0;
	tIBFArbitraryOffsetByPreciseDigraph = 0;
	tLinearUpperRBF = 0;
	tLinearUpperIBF = 0;
	tLinearUpperCSUM = 0;
}

void SchedulabilityAnalysis::reset() {
	// save the results into the vectors
	save_results();
	set_zero();
}

void SchedulabilityAnalysis::output_one_vector(ofstream& fout, string sVec, vector<double> vec) {
	fout<<sVec<<"=[";
	typedef vector<double>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		fout<<*it;
		if (it != --vec.end()) fout<<",";
	}
	fout<<"];"<<endl;
}

void SchedulabilityAnalysis::output_one_vector(ofstream& fout, string sVec, vector<int> vec) {
	fout<<sVec<<"=[";
	typedef vector<int>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		fout<<*it;
		if (it != --vec.end()) fout<<",";
	}
	fout<<"];"<<endl;
}

void SchedulabilityAnalysis::output_vectors(ofstream& fout)  {
	fout<<endl;
	fout<<"Output the informations of Digraphs"<<endl;
	fout<<endl;

	output_one_vector(fout,"tfRatio0",vec_tfRatio0);
	output_one_vector(fout,"tfRatio1",vec_tfRatio1);

	output_one_vector(fout,"nIrreducibleDigraphs",vec_nIrreducibleDigraphs);

	output_one_vector(fout,"tDigraphCalSum",vec_tDigraphCalSum);
	output_one_vector(fout,"tDigraphCalLinearBounds",vec_tDigraphCalLinearBounds);
	output_one_vector(fout,"tDigraphCalTF0",vec_tDigraphCalTF0);
	output_one_vector(fout,"tDigraphCalTF1",vec_tDigraphCalTF1);
	output_one_vector(fout,"tDigraphCalTF2",vec_tDigraphCalTF2);
	output_one_vector(fout,"tDigraphRBF",vec_tDigraphRBF);
	output_one_vector(fout,"tDigraphIBF",vec_tDigraphIBF);

	fout<<endl;
	fout<<"Output the informations of Stateflows"<<endl;
	fout<<endl;

	output_one_vector(fout,"nStateflows",vec_nStateflows);
	output_one_vector(fout,"totalUtil",vec_totalUtilization);
	output_one_vector(fout,"avgDegree",vec_avgDegree);

	output_one_vector(fout,"nExactStaticOffset", vec_nExactStaticOffset);
	output_one_vector(fout,"nExactArbitraryOffset", vec_nExactArbitraryOffset);
	output_one_vector(fout,"nRBFStaticOffset",vec_nRBFStaticOffset);
	output_one_vector(fout,"nRBFArbitraryOffset",vec_nRBFArbitraryOffset);
	output_one_vector(fout,"nRBFArbitraryOffsetBySimpleDigraph",vec_nRBFArbitraryOffsetBySimpleDigraph);
	output_one_vector(fout,"nRBFArbitraryOffsetByPreciseDigraph",vec_nRBFArbitraryOffsetByPreciseDigraph);
	output_one_vector(fout,"nIBFStaticOffset",vec_nIBFStaticOffset);
	output_one_vector(fout,"nIBFArbitraryOffset",vec_nIBFArbitraryOffset);
	output_one_vector(fout,"nIBFArbitraryOffsetBySimpleDigraph",vec_nIBFArbitraryOffsetBySimpleDigraph);
	output_one_vector(fout,"nIBFArbitraryOffsetByPreciseDigraph",vec_nIBFArbitraryOffsetByPreciseDigraph);
	output_one_vector(fout,"nLinearUpperRBF",vec_nLinearUpperRBF);
	output_one_vector(fout,"nLinearUpperIBF",vec_nLinearUpperIBF);
	output_one_vector(fout,"nLinearUpperCSUM", vec_nLinearUpperCSUM);

	output_one_vector(fout,"bpRatio0",vec_bpRatio0);
	output_one_vector(fout,"bpRatio1",vec_bpRatio1);

	output_one_vector(fout,"nIrreducibleStateflows",vec_nIrreducibleStateflows);

	output_one_vector(fout,"tPrepareStateflow",vec_tPrepareStateflow);

	output_one_vector(fout,"tCalCSum",vec_tCalCSum);
	output_one_vector(fout,"tCalLinearBounds",vec_tCalLinearBounds);
	output_one_vector(fout,"tCalTF0",vec_tCalTF0);
	output_one_vector(fout,"tCalTF1",vec_tCalTF1);
	output_one_vector(fout,"tCalTF2",vec_tCalTF2);

	output_one_vector(fout,"tCalNonPeirod",vec_tCalNonPeriod);
	output_one_vector(fout,"tCalPeriod",vec_tCalPeriod);
	output_one_vector(fout,"tCalDiffPeriod",vec_tCalDiffPeriod);

	output_one_vector(fout,"tCalRequestExecutionMatrixWithoutPeriodicityProperty",vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty);
	output_one_vector(fout,"tCalRequestExecutionMatrixWithPeriodicityProperty",vec_tCalRequestExecutionMatrixWithPeriodicityProperty);

	output_one_vector(fout,"tGenerateCriticalPairs",vec_tGenerateCriticalActionPairs);
	output_one_vector(fout,"tGenerateRequestFunctions", vec_tGenerateRequestFunctions);
	output_one_vector(fout,"tGenerateRequestFunctionAbstractTree", vec_tGenerateRequestFunctionAbstractTree);
	
	output_one_vector(fout,"tExactStaticOffset", vec_tExactStaticOffset);
	output_one_vector(fout,"tExactArbitraryOffset", vec_tExactArbitraryOffset);
	output_one_vector(fout,"tRBFStaticOffset",vec_tRBFStaticOffset);
	output_one_vector(fout,"tRBFArbitraryOffset",vec_tRBFArbitraryOffset);
	output_one_vector(fout,"tRBFArbitraryOffsetBySimpleDigraph",vec_tRBFArbitraryOffsetBySimpleDigraph);
	output_one_vector(fout,"tRBFArbitraryOffsetByPreciseDigraph",vec_tRBFArbitraryOffsetByPreciseDigraph);
	output_one_vector(fout,"tIBFStaticOffset",vec_tIBFStaticOffset);
	output_one_vector(fout,"tIBFArbitraryOffset",vec_tIBFArbitraryOffset);
	output_one_vector(fout,"tIBFArbitraryOffsetBySimpleDigraph",vec_tIBFArbitraryOffsetBySimpleDigraph);
	output_one_vector(fout,"tIBFArbitraryOffsetByPreciseDigraph",vec_tIBFArbitraryOffsetByPreciseDigraph);
	output_one_vector(fout,"tLinearUpperRBF",vec_tLinearUpperRBF);
	output_one_vector(fout,"tLinearUpperIBF",vec_tLinearUpperIBF);
	output_one_vector(fout,"tLinearUpperCSUM", vec_tLinearUpperCSUM);
}

void SchedulabilityAnalysis::prepare_all_digraphs(Digraph** digraphs, int n) {
	Timer timer;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();

		timer.start();
		digraph->calculate_csum();
		timer.end();
		tDigraphCalSum += timer.getTime();

		timer.start();
		digraph->calculate_tf0(digraphs,i);
		timer.end();
		tDigraphCalTF0 += timer.getTime();

		timer.start();
		digraph->calculate_linear_upper_bounds();
		timer.end();
		tDigraphCalLinearBounds += timer.getTime();

		timer.start();
		digraph->calculate_tf1(digraphs,i);
		timer.end();
		tDigraphCalTF1 += timer.getTime();

		timer.start();
		digraph->calculate_tf2(digraphs,i);
		timer.end();
		tDigraphCalTF2 += timer.getTime();
		
		tfRatio0 += digraph->tf1/digraph->tf0;
		tfRatio1 += digraph->tf2/digraph->tf0;
	}
}

void SchedulabilityAnalysis::generate_critical_vertices(Digraph** digraphs, int n) {
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		// Determine intra-digraph dominances
		set<Node*> Dominated;

		for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
			Node* node0 = *iter;
			for (vector<Node*>::iterator iter1 = digraph->node_vec.begin(); iter1 != digraph->node_vec.end(); iter1++) {
				Node* node1 = *iter1;
				if (node0 == node1) continue;
				if (Dominated.find(node1) != Dominated.end()) continue;
				if (node0->wcet >= node1->wcet && node0->deadline <= node1->deadline) {
					Dominated.insert(node1);
				}
			}
		}

		for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
			Node* node = *iter;
			if (Dominated.find(node) == Dominated.end()) digraph->cnode_vec.push_back(node);
		}

		// Determine inter-digraph dominances
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			set<Node*> Dominated2;
			vector<Node*> cnode_vec;
			for (vector<Node*>::iterator iiter = digraph->cnode_vec.begin(); iiter != digraph->cnode_vec.end(); iiter++) {
				Node* inode = *iiter;
				for (vector<Node*>::iterator jiter = digraphj->cnode_vec.begin(); jiter != digraphj->cnode_vec.end(); jiter++) {
					Node* jnode = *jiter;
					if (Dominated2.find(jnode) != Dominated2.end()) continue;
					if (inode->wcet >= jnode->wcet && inode->deadline <= jnode->wcet) Dominated2.insert(jnode);
				}
			}

			for (vector<Node*>::iterator iter = digraphj->cnode_vec.begin(); iter != digraphj->cnode_vec.end(); iter++) {
				Node* node = *iter;
				if (Dominated2.find(node) == Dominated2.end()) cnode_vec.push_back(node);
			}

			digraphj->cnode_vec = cnode_vec;
		}
	}
}

bool SchedulabilityAnalysis::rbf_analysis(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();
	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int tf;
		if (choice == 0) tf = ceil(digraph->tf0);
		else if (choice == 1) tf = ceil(digraph->tf1);
		else if (choice == 2) tf = ceil(digraph->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		digraph->prepare_rbf_calculation_without_periodicity(false);
	}

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		int tf;
		if (choice == 0) tf = ceil(digraph->tf0);
		else if (choice == 1) tf = ceil(digraph->tf1);
		else if (choice == 2) tf = ceil(digraph->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf == 0) { schedulable[i] = true; continue; }
		if (tf < digraph->pGCD) tf = digraph->pGCD;

		for (int t=digraph->pGCD; t<=tf; t+=digraph->pGCD) {
			schedulable[i] = false;

			set<int> tprimSet;
			for (int j=0; j<=i; j++) {
				Digraph* hdg = digraphs[j];
				for (int tprim = hdg->pGCD; tprim<=t; tprim+=hdg->pGCD) {
					tprimSet.insert(tprim);
				}
			}
			tprimSet.insert(0);

			for (set<int>::iterator iter = tprimSet.begin(); iter != tprimSet.end(); iter++) {
				int tprim = *iter;
				if (rbf_analysis(digraphs,i,t,tprim)) {
					schedulable[i] = true;
					break;
				}
			}

			if (!schedulable[i]) {
				delete[] schedulable;
				timer.end();
				tDigraphRBF += timer.getTime();
				if (false) cout<<"Failed on Stateflow-"<<i<<" at "<<t<<endl;
				return false;
			}
		}
	}
	
	delete[] schedulable;
	timer.end();
	tDigraphRBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis(Digraph** digraphs, int i, int t, int tprim) {
	double total = 0;
	double dbfi = digraphs[i]->dbf(t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double rbfj = digraphs[j]->rbf(tprim);
		total += rbfj;
		
		if (total > tprim) return false;
	}
	if (total > tprim) return false;
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();
	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int tf;
		if (choice == 0) tf = ceil(digraph->tf0);
		else if (choice == 1) tf = ceil(digraph->tf1);
		else if (choice == 2) tf = ceil(digraph->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		//digraph->prepare_rbf_calculation(false);
		//digraph->tf = tf_max/digraph->gran_digraph->gcd;
		digraph->prepare_ibf_calculation(false);
	}

	

	bool* schedulable = new bool[n];
	for (int i=0; i<n; i++) schedulable[i] = false;

	for (int i=1; i<n; i++) {
		Digraph* digraph = digraphs[i];
		int tf;
		if (choice == 0) tf = ceil(digraph->tf0);
		else if (choice == 1) tf = ceil(digraph->tf1);
		else if (choice == 2) tf = ceil(digraph->tf2);
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		if (tf == 0) { schedulable[i] = true; continue; }
		if (tf < digraph->pGCD) tf = digraph->pGCD;

		for (int t=digraph->pGCD; t<=tf; t+=digraph->pGCD) {
			schedulable[i] = false;

			set<int> tprimSet;
			for (int j=0; j<=i; j++) {
				Digraph* hdg = digraphs[j];
				for (int tprim = hdg->aGCD; tprim<=t; tprim+=hdg->aGCD) {
					tprimSet.insert(tprim);
				}
			}
			tprimSet.insert(0);

			for (set<int>::iterator iter = tprimSet.begin(); iter != tprimSet.end(); iter++) {
				int tprim = *iter;
				if (ibf_analysis(digraphs,i,t,tprim)) {
					schedulable[i] = true;
					break;
				}
			}

			if (!schedulable[i]) {
				delete[] schedulable;
				timer.end();
				tDigraphIBF += timer.getTime();
				return false;
			}
		}
	}
	
	delete[] schedulable;
	timer.end();
	tDigraphIBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis(Digraph** digraphs, int i, int t, int tprim) {
	double total = 0;
	double dbfi = digraphs[i]->dbf(t);
	total += dbfi;

	for (int j=0; j<i; j++) {
		double rbfj = digraphs[j]->ibf(tprim);
		total += rbfj;
		
		if (total > tprim) return false;
	}
	if (total > tprim) return false;
	return true;
}

void SchedulabilityAnalysis::prepare_all_stateflows(Stateflow** sfs, int n, int choice) {
	Timer timer;
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_t_gcd();
		sf->calculate_deadlines();

		//sf->write_graphviz(cout);
		//sf->precise_digraph->write_graphviz(cout);

		//sf->check_irreducible();

		sf->calculate_generialized_period();
		sf->calculate_generialized_defect();
		
		//if (choice == 0) {
			timer.start();
			sf->calculate_csum();
			timer.end();
			tCalCSum += timer.getTime();

			timer.start();
			sf->calculate_tf0(sfs,i);
			timer.end();
			tCalTF0 += timer.getTime();
		//}
		
		//sf->generate_simple_digraph();

		if (choice == 1 || choice == 2) {
			timer.start();
			sf->generate_precise_digraph();
			sf->calculate_linear_upper_bounds(true);
			timer.end();
			tCalLinearBounds += timer.getTime();

			timer.start();
			sf->calculate_tf1(sfs,i);
			timer.end();
			tCalTF1 += timer.getTime();

			timer.start();
			sf->calculate_tf2(sfs,i);
			timer.end();
			tCalTF2 += timer.getTime();

			bpRatio0 += sf->tf1/sf->tf0;
			bpRatio1 += sf->tf2/sf->tf0;
		}
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	// generate all execution request matrices
	calculate_request_execution_matrices_without_periodicity_property(sfs,n,choice);
}

void SchedulabilityAnalysis::reset_calculating_containers(Stateflow** sfs, int n) {
	// reset all the calculating containers
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		
		// release rbf_vec
		for (vector<StartFinishTime*>::iterator iter = sf->rbf_vec.begin(); iter != sf->rbf_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		sf->rbf_vec.clear();

		// release ibf_vec
		for (vector<StartFinishTime*>::iterator iter = sf->ibf_vec.begin(); iter != sf->ibf_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		sf->ibf_vec.clear();
	}
}

void SchedulabilityAnalysis::calculate_request_execution_matrices_without_periodicity_property(Stateflow** sfs, int n, int choice) {
	Timer timer;
	//timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}
	timer.start();
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	timer.end();
	tCalRequestExecutionMatrixWithoutPeriodicityProperty += timer.getTime();
}

void SchedulabilityAnalysis::calculate_request_execution_matrices_with_periodicity_property(Stateflow** sfs, int n, int choice) {
	Timer timer;
	//timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}
	timer.start();
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	timer.end();
	tCalRequestExecutionMatrixWithPeriodicityProperty += timer.getTime();
}

void SchedulabilityAnalysis::generate_critical_action_pair(Stateflow** sfs, int n) {
	Timer timer;
	timer.start();
	int hyperperiod = 1;
	typedef vector<Transition*>::iterator TranVecIter;
	typedef list<Transition*>::iterator TranListIter;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		hyperperiod = Utility::math_lcm(hyperperiod, sf->hyperperiod);
		int maxDeadline = INT_MIN;
		// generate all the action pairs
		for(TranVecIter iter = sf->trans.begin(); iter != sf->trans.end(); iter ++) {
			Transition* tran = *iter;
			// generate all the action pairs for this transition
			for (int stime = 0; stime < hyperperiod; stime += tran->period) {
				State* snk = tran->snk;
				int deadline = INT_MAX;

				if (snk->out.empty()) deadline = stime+tran->period;

				for (TranListIter iter2 = snk->out.begin(); iter2 != snk->out.end(); iter2++) {
					int t = 0;
					while (t<=stime) t+=(*iter2)->period; 
					deadline = min(deadline,t);
				}
				deadline = deadline - stime;
				maxDeadline = max(maxDeadline, deadline);
				ActionPair ap = ActionPair(tran->wcet,deadline,stime,i,tran);
				vector<ActionPair>& vec_AP = sf->mCAP[stime];
				// Determine intra-FSM dominances
				int nSize = vec_AP.size();
				bool* found = new bool[nSize];
				for (int k=0; k<nSize; k++) found[k] = false;
				bool dominated = false;
				bool flag = false;
				int nItem = 0;
				for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
					ActionPair nap = *apIter;
					if (nap>ap) { dominated = true; break;}
					if (ap>nap) { flag = true; found[nItem] = true;}
					nItem ++;
				}

				if (flag) {
					vector<ActionPair> tempAP;
					nItem = 0;
					for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
						ActionPair nap = *apIter;
						if (!found[nItem]) tempAP.push_back(nap); 
						nItem ++;
					}
					tempAP.push_back(ap);
					sf->mCAP[stime] = tempAP;
				}
				else {
					if (!dominated)	
						vec_AP.push_back(ap);
				}

				delete[] found;
			}
		}

		sf->maxDeadline = maxDeadline;

		// Determine inter-FSM dominances
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int,vector<ActionPair>>::iterator mIter;
			for(mIter = sfj->mCAP.begin(); mIter != sfj->mCAP.end(); mIter++) {
				int key = (*mIter).first;
				vector<ActionPair> vec_AP = (*mIter).second;

				if (sf->mCAP.find(key) == sf->mCAP.end()) continue;

				vector<ActionPair> vec_AP2 = sf->mCAP[key];
				int nSize = vec_AP.size();
				bool* found = new bool[nSize];
				for (int k=0; k<nSize; k++) found[k] = false;
				bool flag = false;
				int nItem = 0;

				for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
					ActionPair ap = *apIter;

					for (vector<ActionPair>::iterator apIter2 = vec_AP2.begin(); apIter2 != vec_AP2.end(); apIter2++) {
						ActionPair ap2 = *apIter2;
						if (ap2>ap) { flag = true; found[nItem] = true; break;}
					}
					nItem ++;
				}

				if (flag) {
					vector<ActionPair> tempAP;
					nItem = 0;

					for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
						ActionPair ap = *apIter;

						if (!found[nItem]) tempAP.push_back(ap);
						nItem++;
					}

					sfj->mCAP[key] = tempAP;
				}
				delete[] found;
			}
		}
	}
	timer.end();
	tGenerateCriticalActionPairs += timer.getTime();
}

bool SchedulabilityAnalysis::generate_request_function(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=0; i<n-1; i++) {
		Stateflow* sfi = sfs[i];
		int maxDeadline = INT_MIN;
		int tk = INT_MIN;
		for (int j=i+1; j<n; j++) {
			Stateflow* sfj = sfs[j];
			maxDeadline = max(maxDeadline, sfj->maxDeadline);

			for (map<int,vector<ActionPair>>::iterator mIter = sfj->mCAP.begin(); mIter != sfj->mCAP.end(); mIter++) {
				int first = mIter->first;
				tk = max(tk,first);
			}
		}
		// generate request functions
		bool success = generate_request_function(sfi,i,tk+maxDeadline,output);
		if (!success) {
			timer.end();
			tGenerateRequestFunctions+=timer.getTime();
			return false;
		}
	}

	timer.end();
	tGenerateRequestFunctions+=timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::generate_request_function(Stateflow* sf, int i, int maxDeadline, bool output) {
	if (maxDeadline == 0) maxDeadline = sf->hyperperiod;
	int gcd = sf->gcd;
	int hyperperiod = sf->hyperperiod;
	vector<RequestFunction> vec_rf;

	int t = 0;
	int count = 0;

	while (t <= hyperperiod + maxDeadline) {
		for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter++) {
			int instance = *iter;
			if (instance == hyperperiod) continue;

			t = count*hyperperiod + instance;
			if (t > hyperperiod+maxDeadline) break;

			vector<RequestFunction> temp_vec_rf;
			if (t==0) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
				vec_rf.clear();
				vec_rf = temp_vec_rf;
				temp_vec_rf.clear();
				if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
				continue;
			}

			for (vector<RequestFunction>::iterator rfIter = vec_rf.begin(); rfIter != vec_rf.end(); rfIter++) {
				RequestFunction rf = *rfIter;
				list<Transition*> list_tran = rf.lastTran->snk->out;
				list<Transition*>::iterator tranIter;
				for (tranIter = list_tran.begin(); tranIter != list_tran.end(); tranIter++) {
					RequestFunction temp_rf2 = RequestFunction(rf);
					Transition* tran = *tranIter;
					if (t%tran->period==0) {
						temp_rf2.update_next_action(t,tran);
					} else {
						temp_rf2.update_next_action(t,NULL);
					}
					temp_vec_rf.push_back(temp_rf2);
				}

				if (list_tran.empty()) {
					rf.update_next_action(t, NULL);
					temp_vec_rf.push_back(rf);
				}
			}

			// add other request functions triggered within one hyperperiod but 0
			if (t < hyperperiod) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
			}

			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "=========================" <<endl;
			}
			intra_loop_dominance(temp_vec_rf);
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "+++++++++++++++++++++++++" << endl;
			}
			vec_rf.clear();
			vec_rf = temp_vec_rf;
			temp_vec_rf.clear();
			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			if (vec_rf.size() >= 100) return false;
		}
		count++;
	}

	//output_critical_request_function(vec_rf,cout);
	end_loop_dominance(vec_rf);
	//output_critical_request_function(vec_rf,cout);
	sf->CRF = vec_rf;
	return true;
}

void SchedulabilityAnalysis::intra_loop_dominance(vector<RequestFunction>& vec_rf) {
	// detemine the dominance for rfs
	int nSize = vec_rf.size();
	vector<bool> found = vector<bool>(nSize,false);
	vector<int> count = vector<int>(nSize,0);

	for (int j=0; j<nSize; j++) {
		RequestFunction rf1 = vec_rf.at(j);
		if (found[j]) continue;
		for (int k=0; k<nSize; k++) {
			if (k==j) continue;
			if (found[k]) continue;
			RequestFunction rf2 = vec_rf.at(k);
			if (rf1.start == rf2.start && rf1.lastTran == rf2.lastTran && rf1.isCounted == rf2.isCounted && rf1 > rf2) { 
				found[k] = true;
				count[j]++;
			}
		}
	}

	vector<RequestFunction> temp_vec_rf;
	for (int k=0; k<nSize; k++) {
		if (!found[k]) temp_vec_rf.push_back(vec_rf.at(k));
	}

	if (temp_vec_rf.empty()) {
		int index = 0;
		int temp = INT_MIN;
		for (int k=0; k<nSize; k++) {
			if (count[k] > temp) {
				temp = count[k];
				index = k;
			}
		}
		//cout << vec_rf.size() << "\t" << index <<endl;
		temp_vec_rf.push_back(vec_rf.at(index));
	}

	vec_rf.clear();
	vec_rf = temp_vec_rf;
	temp_vec_rf.clear();

	found.clear();
	count.clear();
}

void SchedulabilityAnalysis::end_loop_dominance(vector<RequestFunction>& vec_rf) {
	// detemine the dominance for rfs
	int nSize = vec_rf.size();
	vector<bool> found = vector<bool>(nSize,false);
	vector<int> count = vector<int>(nSize,0);
	

	for (int j=0; j<nSize; j++) {
		if (found[j]) continue;
		RequestFunction rf1 = vec_rf.at(j);
		for (int k=0; k<nSize; k++) {
			if (k==j) continue;
			if (found[k]) continue;
			RequestFunction rf2 = vec_rf.at(k);
			if (rf1.start == rf2.start && rf1 > rf2) { 
				found[k] = true;
				count[j]++;
			}
		}
	}

	vector<RequestFunction> temp_vec_rf;
	for (int k=0; k<nSize; k++) {
		if (!found[k]) temp_vec_rf.push_back(vec_rf.at(k));
	}

	if (temp_vec_rf.empty()) {
		int index = 0;
		int temp = INT_MIN;
		for (int k=0; k<nSize; k++) {
			if (count[k] > temp) {
				temp = count[k];
				index = k;
			}
		}
		temp_vec_rf.push_back(vec_rf.at(index));
	}

	vec_rf.clear();
	vec_rf = temp_vec_rf;
	temp_vec_rf.clear();

	found.clear();
	count.clear();
}


void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=1; i<n; i++) {
		Stateflow* sfi = sfs[i];
		for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;
			for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
				ActionPair ap = *apIter;
				int tk = 0;
				for (int align = 0; align <=ap.stime; align += sfi->hyperperiod) {
					bool check = true;
					for (int j=0; j<i; j++) {
						Stateflow* sfj = sfs[j];
						if (align % sfj->hyperperiod != 0) {
							check = false;
							break;
						}
					}
					if (check) tk = align;
				}

				tk = ap.stime - tk;

				generate_request_function_abstract_tree(sfs,i,tk,ap.d, output);
			}
		}
	}
	timer.end();
	tGenerateRequestFunctionAbstractTree += timer.getTime();
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow** sfs, int i, int stime, int deadline, bool output) {
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		generate_request_function_abstract_tree(sfj,stime,deadline,output);
	}
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow* sf, int stime, int deadline, bool output) {
	//deadline = ceil(1.0*deadline/sf->gcd)*sf->gcd;

	if (sf->mmARFT.find(stime) != sf->mmARFT.end()) {
		if (sf->mmARFT[stime].find(deadline) != sf->mmARFT[stime].end())
			return;
	}

	AbstractRequestFunctionTree arft = AbstractRequestFunctionTree();
	int gcd = Utility::math_gcd(sf->gcd,deadline);
	arft.generate(sf->CRF,stime,deadline,gcd,output);

	sf->mmARFT[stime][deadline] = arft;
}

bool SchedulabilityAnalysis::isAbstract(vector<RFNode*> vec_rfnode) {
	for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
		RFNode* rfnode = *iter;
		if (rfnode->isAbstract) return true;
	}
	return false;
}

bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int n, int choice, bool output, ostream& out) {
	Timer timer;
	timer.start();

	if (output) cout << "Generating request functions" << endl;
	bool success = SchedulabilityAnalysis::generate_request_function(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_critical_request_function(sfs,n,out);

	if (output) cout << "Generating request function abstract tree" <<endl; 
	SchedulabilityAnalysis::generate_request_function_abstract_tree(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_request_function_abstract_tree(sfs,n,out);

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		for (map<int, vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			bool schedulable = exact_sched_analysis(sfs,i,stime,vec_ap);
			if (!schedulable) {
				timer.end();
				tExactStaticOffset += timer.getTime();
				return false;
			}
		}
	}
	timer.end();
	tExactStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int i, int s, vector<ActionPair> vec_ap) {
	for (vector<ActionPair>::iterator iter = vec_ap.begin(); iter != vec_ap.end(); iter++) {
		bool schedulable = exact_sched_analysis(sfs, i, s, *iter);
		if (!schedulable) return false;
	}
	return true;
}


bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int i, int s, ActionPair ap) {
	if (i==0) {
		if (ap.e <= ap.d) return true;
		else return false;
	} else {
		int gcd = 0;
		for (int j=0; j<=i; j++) gcd = Utility::math_gcd(gcd, sfs[j]->gcd);

		vector<AbstractRequestFunctionTree> vec_arft;

		int tk = 0;
		for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
			bool check = true;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				if (align % sfj->hyperperiod != 0) {
					check = false;
					break;
				}
			}
			if (check) tk = align;
		}
		tk = ap.stime - tk;
		
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int, map<int,AbstractRequestFunctionTree>> mmARFT = sfj->mmARFT;
			if (mmARFT.find(tk) == mmARFT.end()) {
				cerr << "Not generate " << tk << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}

			if (mmARFT[tk].find(ap.d) == mmARFT[tk].end()) {
				cerr << "Not generate " << tk << "," << ap.d << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}
			AbstractRequestFunctionTree arft = mmARFT[tk][ap.d];
			vec_arft.push_back(arft);
		}

		return exact_sched_analysis(vec_arft,ap.e,tk,ap.d,gcd);
	}
}


bool SchedulabilityAnalysis::exact_sched_analysis(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int stime, int deadline, int gcd) {
	vector<vector<RFNode*>> Omega;
	vector<RFNode*> vec_rf;
	for (vector<AbstractRequestFunctionTree>::iterator iter = vec_arft.begin(); iter != vec_arft.end(); iter++) {
		vec_rf.push_back((*iter).root);
	}
	Omega.push_back(vec_rf);

	return exact_sched_analysis(Omega, wcet, stime, deadline, gcd);
}

bool SchedulabilityAnalysis::exact_sched_analysis(vector<vector<RFNode*>> Omega, int wcet, int stime, int deadline, int gcd) {
	for (vector<vector<RFNode*>>::iterator iter = Omega.begin(); iter != Omega.end(); iter++) {
		vector<RFNode*> vec_rfnode = *iter;
		bool schedulable = exact_sched_analysis(vec_rfnode, wcet, stime, deadline, gcd);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_sched_analysis(vector<RFNode*> vec_rfnode, int wcet, int stime, int deadline, int gcd) {
	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);

	// check the schedulability condition
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		if (rt == INT_MAX) break;
		int start = *iter;
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double sum = 0;
			if (tprim >= stime) sum += wcet;
			for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
				RFNode* rfnode = *iter;
				sum += rfnode->get_mrf(start,tprim);
			}
			if (sum <= tprim-start) {
				rt = max(rt,sum+start-stime);
				break;
			}
			if (tprim == stime+deadline) rt=INT_MAX;
		}
	}

	if (rt <= deadline) return true;

	// we need to consider the unschedulable situation
	if (isAbstract(vec_rfnode)) {
		vector<RFNode*> vec_rfnode1;
		vector<RFNode*> vec_rfnode2;

		bool found = false;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			if (rfnode->isAbstract && !found) {
				found = true;
				vec_rfnode1.push_back(rfnode->lnode);
				vec_rfnode2.push_back(rfnode->rnode);
				continue;
			}
			vec_rfnode1.push_back(rfnode);
			vec_rfnode2.push_back(rfnode);
		}

		bool schedulable1 = exact_sched_analysis(vec_rfnode1,wcet, stime, deadline, gcd);
		bool schedulable2 = exact_sched_analysis(vec_rfnode2,wcet, stime, deadline, gcd);
		
		return schedulable1 && schedulable2;

	} else 
		return false;
}

int SchedulabilityAnalysis::calculate_remaining_workload(vector<RFNode*> vec_rfnode, int stime, int gcd) {
	int remain = 0;
	int interval = 0;
	int prev = 0;
	for (int t=gcd; t<=stime; t+=gcd) {
		interval += gcd;
		int curr = 0;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			if (rfnode->stime != stime) {
				cerr << "calculating remaining workload error" << rfnode->stime << "\t" << stime <<endl;
				exit(EXIT_FAILURE);
			}
			curr += rfnode->get_mrf(0,t);
		}
		remain += curr-prev;
		prev = curr;
		if (remain <= interval) {
			remain = 0;
			interval = 0;
		}
	}

	return remain;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	// release critical request functions
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		// release abstract request function trees
		for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mmIter = sf->mmARFT.begin(); mmIter != sf->mmARFT.end(); mmIter++) {
			map<int,AbstractRequestFunctionTree> temp = mmIter->second;
			for (map<int,AbstractRequestFunctionTree>::iterator mIter = temp.begin(); mIter != temp.end(); mIter++) {
				AbstractRequestFunctionTree& arft = mIter->second;
				arft.root = NULL;
				while(!arft.cur_queue.empty()) arft.cur_queue.pop();
				for (vector<RFNode*>::iterator iter = arft.record.begin(); iter != arft.record.end(); iter++) {
					delete *iter;
				}
				arft.record.clear();
			}
		}

		sf->CRF.clear();
		sf->mmARFT.clear();
	}
	
	Timer timer;
	timer.start();

	bool output = false;

	if (output) cout << "Generating request functions" << endl;
	bool success = SchedulabilityAnalysis::generate_request_function(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_critical_request_function(sfs,n,cout);

	if (output) cout << "Generating request function abstract tree" <<endl; 
	SchedulabilityAnalysis::generate_request_function_abstract_tree(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_request_function_abstract_tree(sfs,n,cout);

	// request functions and abstract rf trees have been generated in function exact_sched_analysis
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		for (map<int, vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			bool schedulable = exact_analysis_arbitrary_offset(sfs,i,vec_ap);
			if (!schedulable) {
				timer.end();
				tExactArbitraryOffset += timer.getTime();
				return false;
			}
		}
	}
	timer.end();
	tExactArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int i, vector<ActionPair> vec_ap) {
	for (vector<ActionPair>::iterator iter = vec_ap.begin(); iter != vec_ap.end(); iter++) {
		bool schedulable = exact_analysis_arbitrary_offset(sfs, i, *iter);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int i, ActionPair ap) {
	if (i==0) {
		if (ap.e <= ap.d) return true;
		else return false;
	} else {
		int gcd = 0;
		for (int j=0; j<=i; j++) gcd = Utility::math_gcd(gcd, sfs[j]->gcd);

		vector<AbstractRequestFunctionTree> vec_arft;
		int tk =ap.stime;

		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int, map<int,AbstractRequestFunctionTree>> mmARFT = sfj->mmARFT;
			if (mmARFT.find(tk) == mmARFT.end()) {
				cerr << "Not generate " << tk << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}

			if (mmARFT[tk].find(ap.d) == mmARFT[tk].end()) {
				cerr << "Not generate " << tk << "," << ap.d << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}
			AbstractRequestFunctionTree arft = mmARFT[tk][ap.d];
			vec_arft.push_back(arft);
		}

		return exact_analysis_arbitrary_offset(vec_arft, ap.e,ap.d,gcd);
	}
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int deadline, int gcd) {
	vector<vector<RFNode*>> Omega;
	vector<RFNode*> vec_rf;
	for (vector<AbstractRequestFunctionTree>::iterator iter = vec_arft.begin(); iter != vec_arft.end(); iter++)
		vec_rf.push_back((*iter).root);
	Omega.push_back(vec_rf);
	
	return exact_analysis_arbitrary_offset(Omega, wcet, deadline, gcd);
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<vector<RFNode*>> Omega, int wcet, int deadline, int gcd) {
	for (vector<vector<RFNode*>>::iterator iter = Omega.begin(); iter != Omega.end(); iter++) {
		vector<RFNode*> vec_rfnode = *iter;
		bool schedulable = exact_analysis_arbitrary_offset(vec_rfnode, wcet, deadline, gcd);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<RFNode*> vec_rfnode, int wcet, int deadline, int gcd) {
	// check the schedulability condition
	double rt = 0;
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			sum += rfnode->get_mrf(t);
		}
		if (sum <= t) return true;
	}
	return false;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	/*
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) return false;
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) return false;
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) return false;
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	*/

	for (int i=0; i<n; i++) {
		bool schedulable = rbf_analysis_static_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tRBFStaticOffset += timer.getTime();	
			return false;
		}
	}
	
	timer.end();
	tRBFStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;

			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			bool schedulable = rbf_analysis_static_offset_index(sfs,i,ap.e,tk,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale by rbf analysis with static offset on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;	
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);
	
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				double rbfsf = sfj->get_rbf(start,tprim);
				temp += rbfsf;
			}
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				break;
			}
			if (tprim == stime+deadline) return false;
		}
	}

	if (rt <= deadline) return true;
}

double SchedulabilityAnalysis::calculate_remaining_rbf_workload(Stateflow** sfs, int i, int stime, int gcd) {
	double remain = 0;
	int interval = 0;
	double prev = 0;
	for (int t=gcd; t<=stime; t+=gcd) {
		interval += gcd;
		double curr = 0;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			curr += sfj->get_rbf(0,t);
		}
		remain += curr-prev;
		prev = curr;
		if (remain <= interval) {
			remain = 0;
			interval = 0;
		}
	}

	return remain;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	/*
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) return false;
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) return false;
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) return false;
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	*/

	for (int i=0; i<n; i++) {
		bool schedulable = ibf_analysis_static_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tIBFStaticOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tIBFStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;

			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			bool schedulable = ibf_analysis_static_offset_index(sfs,i,ap.e,tk,ap.d);
			if (!schedulable) return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	int t_gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
		t_gcd = Utility::math_gcd(t_gcd,sfj->t_gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);
	
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		bool found = false;
		int tprim = start+gcd;
		while (tprim <= stime+deadline) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				double ibfsf = sfj->get_ibf(start,tprim);
				temp += ibfsf;
			}
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				found = true;
				break;
			}
			int total = ceil(temp)+start;
			if (total == tprim) tprim += t_gcd;
			tprim = total;
		}
		if (!found) return false;
	}

	if (rt <= deadline) return true;
	
}

double SchedulabilityAnalysis::calculate_remaining_ibf_workload(Stateflow** sfs, int i, int stime, int gcd) {
	return calculate_remaining_rbf_workload(sfs,i,stime,gcd);
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	/*
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) return false;
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) return false;
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) return false;
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}
	*/

	for (int i=0; i<n; i++) {
		bool schedulable = rbf_analysis_arbitrary_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tRBFArbitraryOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tRBFArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			bool schedulable = rbf_analysis_arbitrary_offset_index(sfs,i,ap.e,ap.d);
			if (!schedulable) {
				cout << "Unschedulable by RBF analysis with arbitrary offset for Stateflow " << i 
					<< "\t Action Pair=>" << ap.toString() <<endl;	
				return false;
			}
		}
	}
	
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			sum += sfj->get_rbf(t);
		}
		if (sum <= t) return true;
	}

	return false;
}



bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
	}
	//timer.end();
	//tRBFArbitraryOffsetBySimpleDigraph += timer.getTime();

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,n);

	
	// set the gcds
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->simple_digraph->pGCD = sf->gcd;
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		//digraph->prepare_rbf_calculation_without_periodicity(false);
		digraph->prepare_rbf_calculation(false);
	}
	//timer.start();
	
	//digraphs[6]->write_graphviz(cout);

	bool schedulable = false;

	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(digraphs,i);
		if (!schedulable) break;
	}

	delete[] digraphs;

	timer.end();
	tRBFArbitraryOffsetBySimpleDigraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = rbf_analysis_arbitrary_offset_by_simple_digraphs(digraphs, i, node->wcet, node->deadline);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->pGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			sum += digraphj->rbf(t);
		}
		if (sum <= t) return true;
	}
	return false;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,n);
	bool ret = SchedulabilityAnalysis::rbf_analysis(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	tRBFArbitraryOffsetByPreciseDigraph+= timer.getTime();
	return ret;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	/*
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) return false;
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) return false;
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) return false;
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+1);
	}
	*/

	for (int i=0; i<n; i++) {
		Stateflow* sfi = sfs[i];
		sfi->calculate_deadlines();
		bool schedulable = ibf_analysis_arbitrary_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tIBFArbitraryOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tIBFArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (vector<Transition*>::iterator iter = sfi->trans.begin(); iter != sfi->trans.end(); iter++) {
		Transition* t = *iter;
		bool schedulable = ibf_analysis_arbitrary_offset_index(sfs,i,t->wcet,t->deadline);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->t_gcd);
	}
	
	// check the schedulability condition
	int t = 0;
	while (t<=deadline) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			sum += sfj->get_ibf(t);
		}
		if (sum <= t) return true;
		int total = ceil(sum);
		if (total==t) t+=gcd;
		t = total;
	}

	return false;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,n);

	// set the gcds
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->simple_digraph->pGCD = sf->gcd;
	}

	bool ret = SchedulabilityAnalysis::ibf_analysis(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	tIBFArbitraryOffsetBySimpleDigraph += timer.getTime();
	return ret;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,n);
	bool ret = SchedulabilityAnalysis::ibf_analysis(digraphs,n,choice);

	delete[] digraphs;
	timer.end();
	tIBFArbitraryOffsetByPreciseDigraph += timer.getTime();
	return ret;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	for (int i=0; i<n; i++) {
		bool schedulable = lu_rbf_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperRBF += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperRBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			/*
			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			*/
			bool schedulable = lu_rbf_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	// generate the start times for all the busy periods
	/*
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);
	*/

	double tUtil = 0;
	double tCrbf = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCrbf += sfj->crbf;
		if (tUtil >= 1) return false;
	}

	/*
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			temp += tUtil*(tprim-start)+tCrbf;
			if (start+temp>stime+deadline) return false;
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				break;
			}
			if (tprim == stime+deadline) return false;
		}
	}
	*/
	/*
	for (int tprim = gcd; tprim <= deadline; tprim+=gcd) {
		double temp = 0;
		if (tprim >= stime) temp += wcet;
		temp += tUtil*tprim+tCrbf;
		if (temp>deadline) return false;
		if (temp <= tprim) {
			return true;
		}
	}
	*/

	double rt = (tCrbf+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	for (int i=0; i<n; i++) {
		bool schedulable = lu_ibf_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperIBF += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperIBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			/*
			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			*/
			bool schedulable = lu_ibf_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	/*
	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);

	double tUtil = 0;
	double tCibf = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCibf += sfj->cibf;
		if (tUtil > 1) return false;
	}

	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			temp += tUtil*(tprim-start)+tCibf;
			if (start+temp>stime+deadline) return false;
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				break;
			}
			if (tprim == stime+deadline) return false;
		}
	}

	if (rt <= deadline) return true;
	return false;
	*/

	double tUtil = 0;
	double tCibf = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCibf += sfj->cibf;
		if (tUtil >= 1) return false;
	}

	/*
	for (int tprim = gcd; tprim <= deadline; tprim+=gcd) {
		double temp = 0;
		if (tprim >= stime) temp += wcet;
		temp += tUtil*tprim+tCibf;
		if (temp>deadline) return false;
		if (temp <= tprim) {
			return true;
		}
	}
	return false;
	*/

	double rt = (tCibf+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	for (int i=0; i<n; i++) {
		bool schedulable = lu_csum_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperCSUM += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperCSUM += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			/*
			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			*/
			bool schedulable = lu_csum_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}
	/*
	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);

	double tUtil = 0;
	double tCSum = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCSum += 2.0*sfj->csum;
		if (tUtil > 1) return false;
	}

	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			temp += tUtil*(tprim-start)+tCSum;
			if (start+temp>stime+deadline) return false;
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				break;
			}
			if (tprim == stime+deadline) return false;
		}
	}

	if (rt <= deadline) return true;
	return false;
	*/

	double tUtil = 0;
	double tCSum = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCSum += 2.0*sfj->csum;
		if (tUtil >= 1) return false;
	}
	/*
	for (int tprim = gcd; tprim <= deadline; tprim+=gcd) {
		double temp = 0;
		if (tprim >= stime) temp += wcet;
		temp += tUtil*tprim+tCSum;
		if (temp>deadline) return false;
		if (temp <= tprim) {
			return true;
		}
	}
	return false;
	*/

	double rt = (tCSum+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

void SchedulabilityAnalysis::output_critical_action_pair(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n; i++) {
		out << "Stateflow "<<i<< " "": Show critical action pairs:"<<endl;
		Stateflow* sf = sfs[i];

		for (map<int,vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			out << "\t" << "Trigger time = "<<stime<<endl;
			for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
				ActionPair ap = *apIter;
				out << "\t\t" << "<" << ap.priority << "," << ap.stime << "," << ap.e << "," << ap.d << ">" <<endl; 
			}
		}
		out <<endl;
	}
}

void SchedulabilityAnalysis::output_critical_request_function(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n; i++) {
		out << "Stateflow " <<i<<": Show critical request functions:"<<endl;
		Stateflow* sf = sfs[i];

		output_critical_request_function(sf,out);
	}
}

void SchedulabilityAnalysis::output_critical_request_function(Stateflow* sf, ostream& out) {
	out << "Number=" << sf->CRF.size() <<endl;
	output_critical_request_function(sf->CRF,out);
	out <<endl;
}

void SchedulabilityAnalysis::output_critical_request_function(vector<RequestFunction> vec_rf, ostream& out) {
	for (vector<RequestFunction>::iterator rfIter = vec_rf.begin(); rfIter != vec_rf.end(); rfIter++) {
		RequestFunction rf = *rfIter;
		rf.output(out);
	}

}

void SchedulabilityAnalysis::output_request_function_abstract_tree(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n-1; i++) {
		out << "Stateflow " <<i<<": Show request functions abstract tree:"<<endl;
		Stateflow* sf = sfs[i];

		output_request_function_abstract_tree(sf,out);
	}
}

void SchedulabilityAnalysis::output_request_function_abstract_tree(Stateflow* sf, ostream& out) {
	for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mIter = sf->mmARFT.begin(); mIter != sf->mmARFT.end(); mIter++) {
		int stime = mIter->first;
		map<int,AbstractRequestFunctionTree> mARFT = mIter->second;
		for (map<int,AbstractRequestFunctionTree>::iterator mIter2 = mARFT.begin(); mIter2 != mARFT.end(); mIter2++) {
			int deadline = mIter2->first;
			AbstractRequestFunctionTree arft = mIter2->second;
			out << "Within [" << stime << "," << deadline <<")" <<endl;
			output_request_function_abstract_tree(arft,out);
			out << endl;
		}
	}
}

void SchedulabilityAnalysis::output_request_function_abstract_tree(AbstractRequestFunctionTree arft, ostream& out) {
	arft.root->output(out);
}