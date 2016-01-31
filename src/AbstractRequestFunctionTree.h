#ifndef ABSTRACTREQUESTFUNCTIONTREE_H_
#define ABSTRACTREQUESTFUNCTIONTREE_H_

#include <queue>
#include "RequestFunction.h"

class RFNode {
public:
	RequestFunction rf;
	vector<RequestFunction> vec_rf;
	int gcd;
	int stime;
	int deadline;
	bool isAbstract;
	RFNode* lnode;
	RFNode* rnode;
	map<int,map<int,int>> mrf; // mrf[s,f)
	map<int,int> mrf2; // mrf(t)

	RFNode(const RFNode& rfnode) {
		rf = rfnode.rf;
		gcd = rfnode.gcd;
		stime = rfnode.stime;
		deadline = rfnode.deadline;
		isAbstract = rfnode.isAbstract;
		lnode = rfnode.lnode;
		rnode = rfnode.rnode;
		mrf = rfnode.mrf;
		vec_rf = rfnode.vec_rf;
	}

	RFNode(RequestFunction _rf, int _stime, int _deadline, int _gcd) {
		rf = _rf;
		vec_rf.push_back(rf);

		stime = _stime;
		deadline = _deadline;
		gcd = _gcd;

		isAbstract = false;
		lnode = NULL;
		rnode = NULL;
	}

	RFNode(int _stime, int _deadline, int _gcd, RFNode* _lnode, RFNode* _rnode) {
		for (vector<RequestFunction>::iterator iter = _lnode->vec_rf.begin(); iter != _lnode->vec_rf.end(); iter++) {
			vec_rf.push_back(*iter);
		}

		for (vector<RequestFunction>::iterator iter = _rnode->vec_rf.begin(); iter != _rnode->vec_rf.end(); iter++) {
			vec_rf.push_back(*iter);
		}

		stime = _stime;
		deadline = _deadline;
		gcd = _gcd;

		isAbstract = true;
		lnode = _lnode;
		rnode = _rnode;

		int start = ceil(1.0*stime/gcd)*gcd;

		if (_lnode->stime != stime || _rnode->stime != stime) {
			cerr << "RFNode " << _lnode->stime <<"\t"<<_rnode->stime<<"\t"<<stime<<endl;
			exit(EXIT_FAILURE);
		}
	}

	~RFNode() {
		lnode = NULL;
		rnode = NULL;
		mrf.clear();
	}

	int get_mrf(int s, int f) {
		if ( s > f || s<0) {
			cerr << "Error interval stime="<<stime<<"\t"<<"["<<s<<","<<f<<")"<<endl;
			exit(EXIT_FAILURE);
		}

		s = ceil(1.0*s/gcd)*gcd;
		f = ceil(1.0*f/gcd)*gcd;

		if (s==f) return 0;

		if (mrf.find(s) != mrf.end()) {
			if (mrf[s].find(f) != mrf[s].end())
				return mrf[s][f];
		}

		int value = INT_MIN;
		for (vector<RequestFunction>::iterator iter = vec_rf.begin(); iter != vec_rf.end(); iter++) {
			RequestFunction rf = *iter;
			int rfsf = rf.request_function(s,f);
			value = max(value,rfsf);
		}
		mrf[s][f] = value;
		return value;
	}

	int get_mrf(int t) {
		t = ceil(1.0*t/gcd)*gcd;

		if (mrf2.find(t) != mrf2.end())
			return mrf2[t];

		int value = INT_MIN;
		for (vector<RequestFunction>::iterator iter = vec_rf.begin(); iter != vec_rf.end(); iter++) {
			RequestFunction rf = *iter;
			int rft = rf.request_function(t);
			value = max(value,rft);
		}
		
		mrf2[t] = value;
		return value;
	}

	void output(ostream& out) {
		int s = ceil(1.0*stime/gcd)*gcd;
		int f = ceil(1.0*deadline/gcd)*gcd;

		for (int t = s; t<=s+f; t+=gcd) {
			out<< "mrf[" << s << "," << t << ")=" << get_mrf(s,t) <<endl;
		}
	}
};

class AbstractRequestFunctionTree {
public:
	RFNode* root;
	vector<RFNode*> record;
	queue<RFNode*> cur_queue;

	AbstractRequestFunctionTree () {}

	AbstractRequestFunctionTree (RFNode* _root) {
		root = _root;
		cur_queue.push(root);
	}

	/*
	~AbstractRequestFunctionTree() {
		root = NULL;
		while(!cur_queue.empty()) cur_queue.pop();
		for (vector<RFNode*>::iterator iter = record.begin(); iter != record.end(); iter++) {
			delete *iter;
		}
		record.clear();
	}
	*/
	void refine() {
		RFNode* front = cur_queue.front();
		cur_queue.pop();
		cur_queue.push(front->lnode);
		cur_queue.push(front->rnode);
	}

	void reset() {
		while (!cur_queue.empty()) cur_queue.pop();
		cur_queue.push(root);
	}

	// Similarity metric
	double calculate_dist(RFNode* rfnode1, RFNode* rfnode2, int stime, int deadline, int gcd) {
		int d = ceil(1.0*deadline/gcd);
		double alpha = pow(0.1,1.0/d);
		double sum = 0.0;

		if (rfnode1->stime != stime || rfnode2->stime != stime) {
			cerr << "Calculating distance error " << rfnode1->stime <<"\t"<<rfnode2->stime<<"\t"<<stime<<endl;
			exit(EXIT_FAILURE);
		}

		// front part
		int diff = abs(rfnode1->get_mrf(0,stime)-rfnode2->get_mrf(0,stime));
		sum += pow(alpha,0)*diff;

		// back part
		int i=0;
		for (int t=ceil(stime/gcd)*gcd; t<=ceil((stime+deadline)/gcd)*gcd; t+=gcd) {
			int value1 = rfnode1->get_mrf(stime,t);
			int value2 = rfnode2->get_mrf(stime,t);
			sum += pow(alpha, i)*abs(value1-value2);
		}
		return sum;
	}

	void iterative_generate(vector<RFNode*> vec_rfnode, int stime, int deadline, int gcd, bool output) {
		int nSize = vec_rfnode.size();
		if (output) cout << "nSIze=" << nSize <<endl;
		
		if (nSize == 1) {
			root = vec_rfnode.at(0);
			cur_queue.push(root);
			if (output) cout << record.size() <<endl;
			return;
		}

		double** dist = new double*[nSize];
		for (int i=0; i<nSize; i++) dist[i] = new double[nSize];

		for (int i=0; i<nSize; i++) {
			RFNode* rfnode1 = vec_rfnode.at(i);
			for (int j=0; j<nSize; j++) {
				if (i==j) { 
					dist[i][j] = INT_MAX;
					continue;
				}
				RFNode* rfnode2 = vec_rfnode.at(j);
				dist[i][j] = dist[j][i] = calculate_dist(rfnode1,rfnode2,stime,deadline,gcd);
			}
		}

		queue<int> iQueue;
		set<int> iSet;
		while (iSet.size() != nSize) {
			if(output) cout << "iSet.size=" <<iSet.size() <<endl;
			// find i,j such that the minimal similarity metric
			double minimum = INT_MAX;
			bool found = false;
			int tempX = 0;
			int tempY = 0;
			for (int i=0; i<nSize; i++) {
				if (iSet.find(i) != iSet.end()) continue;
				for (int j=0; j<nSize; j++) {
					if (iSet.find(j) != iSet.end()) continue;
					if (i > j) {
						if (minimum > dist[i][j]) {
							found = true;
							minimum = dist[i][j];
							tempX = i;
							tempY = j;
						}
					}
				}
			}
			if (found) {
				if (tempX == tempY) {
					cerr << "tempX=tempY=" << tempX <<endl;
					exit(EXIT_FAILURE);
				}
				iQueue.push(tempX);
				iQueue.push(tempY);
				iSet.insert(tempX);
				iSet.insert(tempY);

				dist[tempX][tempY] = dist[tempY][tempX] = INT_MAX;
				
			} else {
				if (iSet.size() != nSize-1) {
					cerr << "iSet.size=" << iSet.size() <<"\tnSize="<<nSize<<endl;
					exit(EXIT_FAILURE);
				}
				for (int i=0; i<nSize; i++) {
					if (iSet.find(i) == iSet.end()) {
						iQueue.push(i);
						iSet.insert(i);
						break;
					}
				}
			}
		}

		vector<RFNode*> vec_rfnode2;
		while (!iQueue.empty()) {
			if (iQueue.size() == 1) {
				int index = iQueue.front();
				iQueue.pop();
				vec_rfnode2.push_back(vec_rfnode.at(index));
				break;
			}

			// the size of iQueue should not be less than 2
			int index0 = iQueue.front();
			iQueue.pop();
			int index1 = iQueue.front();
			iQueue.pop();
			RFNode* combine = new RFNode(stime, deadline, gcd,vec_rfnode.at(index0),vec_rfnode.at(index1));
			record.push_back(combine);
			vec_rfnode2.push_back(combine);
		}

		for (int i=0; i<nSize; i++) delete[] dist[i];
		delete[] dist;

		iterative_generate(vec_rfnode2, stime, deadline, gcd,output);
	}

	void generate(vector<RequestFunction> vec_RF, int stime, int deadline, int gcd, bool output) {
		vector<RFNode*> vec_rfnode;
		for (vector<RequestFunction>::iterator iter = vec_RF.begin(); iter != vec_RF.end(); iter++) {
			RFNode* rfnode = new RFNode(*iter,stime,deadline,gcd);
			record.push_back(rfnode);
			vec_rfnode.push_back(rfnode);
		}

		iterative_generate(vec_rfnode, stime, deadline, gcd, output);
	}

	void output(ostream& out) {
		
	}
};

#endif