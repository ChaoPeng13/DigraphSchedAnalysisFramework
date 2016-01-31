#include "FileWriter.h"

void FileWriter::DotFileWriter(Digraph** digraphs, int num, const char* fname) {
	ofstream fout(fname, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		digraph->write_graphviz(fout);
		fout<<endl;
	}
	fout.close();
}

void FileWriter::DotFileWriter(Stateflow* sf, const char* fname) {
	ofstream fout(fname, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	sf->write_graphviz(fout);
	fout<<endl;
	fout.close();
}

void FileWriter::DotFileWriter(Stateflow** sfs, int num, const char* fname) {
	ofstream fout(fname, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		sf->write_graphviz(fout);
		fout<<endl;
	}
	fout.close();
}

void FileWriter::DotFileWriter(Stateflow* sf, int range, string prefix) {
	string file = "TestResult\\"+prefix+".m";
	const char *filename = file.c_str();
	ofstream fout(filename, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<filename<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	fout<<"function "<<prefix<<"()"<<endl;
	fout<<endl;

	// Note that Stateflow sf should be prepared before.
	// write rbf[s,f)
	fout<<"%=================rbf[s,f)======================"<<endl;
	for (set<int>::iterator sIter = sf->rbf_time_instances.begin(); sIter != sf->rbf_time_instances.end(); sIter++) {
		int s = *sIter;
		// write x-direction
		fout<<"x"<<s<<"="<<s<<":1:"<<s+range<<";"<<endl;
		fout<<endl;

		fout<<"rbf"<<s<<"=[";
		for (int f=0; f<=range; f++) {
			double value = sf->get_rbf(s,s+f);
			if (value == NEG_INFINITY) value = 0;
			fout<<value;
			if (f!=range) fout<<",";
		}
		fout<<"];"<<endl;
		fout<<endl;
	}

	// write ibf[s,f)
	fout<<"%=================ibf[s,f)======================"<<endl;
	for (set<int>::iterator sIter = sf->rbf_time_instances.begin(); sIter != sf->rbf_time_instances.end(); sIter++) {
		int s = *sIter;
		fout<<"ibf"<<s<<"=[";
		for (int f=0; f<=range; f++) {
			double value = sf->get_ibf(s,s+f);
			if (value == NEG_INFINITY) value = 0;
			fout<<value;
			if (f!=range) fout<<",";
		}
		fout<<"];"<<endl;
		fout<<endl;
	}

	// write dbf[s,f]
	fout<<"%=================dbf[s,f]======================"<<endl;
	for (set<int>::iterator sIter = sf->rbf_time_instances.begin(); sIter != sf->rbf_time_instances.end(); sIter++) {
		int s = *sIter;
		fout<<"dbf"<<s<<"=[";
		for (int f=0; f<=range; f++) {
			double value = sf->get_dbf(s,s+f);
			if (value == NEG_INFINITY) value = 0;
			fout<<value;
			if (f!=range) fout<<",";
		}
		fout<<"];"<<endl;
		fout<<endl;
	}

	// write rbf(t)
	fout<<"%=================rbf(t)======================"<<endl;
	fout<<"rbf=[";
	for (int t=0; t <= range; t++) {
		fout<<sf->get_rbf(t);
		if (t!=range) fout<<",";
	}
	fout<<"];"<<endl;
	fout<<endl;

	// write ibf(t)
	fout<<"%=================ibf(t)======================"<<endl;
	fout<<"ibf=[";
	for (int t=0; t <= range; t++) {
		fout<<sf->get_ibf(t);
		if (t!=range) fout<<",";
	}
	fout<<"];"<<endl;
	fout<<endl;

	// write dbf(t)
	fout<<"%=================dbf(t)======================"<<endl;
	fout<<"dbf=[";
	for (int t=0; t <= range; t++) {
		fout<<sf->get_dbf(t);
		if (t!=range) fout<<",";
	}
	fout<<"];"<<endl;
	fout<<endl;

	// plot these lines
	fout<<"%=================Draw Graphs======================"<<endl;
	fout<<"plot"<<endl;

	fout.close();

}

void FileWriter::WriteRbfAndDbf(Stateflow** sfs, int num, ostream& out) {
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		out<<"Stateflow-"<<i<<endl;

		out<<"%=================rbf[s,f)======================"<<endl;
		for (set<int>::iterator sIter = sf->rbf_time_instances.begin(); sIter != sf->rbf_time_instances.end(); sIter++) {
			int s = *sIter;
			// write x-direction
			out<<"x"<<s<<"="<<s<<":"<<sf->gcd<<":"<<s+sf->hyperperiod<<";"<<endl;
			out<<endl;

			out<<"rbf"<<s<<"=[";
			for (int f=0; f<=sf->hyperperiod; f+=sf->gcd) {
				double value = sf->get_rbf(s,s+f);
				if (value == NEG_INFINITY) value = 0;
				out<<value;
				if (f!=sf->hyperperiod) out<<",";
			}
			out<<"];"<<endl;
			out<<endl;
		}

		// write dbf[s,f]
		out<<"%=================dbf[s,f]======================"<<endl;
		for (set<int>::iterator sIter = sf->rbf_time_instances.begin(); sIter != sf->rbf_time_instances.end(); sIter++) {
			int s = *sIter;
			out<<"dbf"<<s<<"=[";
			for (int f=0; f<=sf->hyperperiod; f+=sf->gcd) {
				double value = sf->get_dbf(s,s+f);
				if (value == NEG_INFINITY) value = 0;
				out<<value;
				if (f!=sf->hyperperiod) out<<",";
			}
			out<<"];"<<endl;
			out<<endl;
		}

		// write rbf(t)
		out<<"%=================rbf(t)======================"<<endl;
		out<<"rbf=[";
		for (int t=0; t <= sf->hyperperiod; t+=sf->gcd) {
			out<<sf->get_rbf(t);
			if (t!=sf->hyperperiod) out<<",";
		}
		out<<"];"<<endl;
		out<<endl;

		// write dbf(t)
		out<<"%=================dbf(t)======================"<<endl;
		out<<"dbf=[";
		for (int t=0; t <= sf->hyperperiod; t+=sf->gcd) {
			out<<sf->get_dbf(t);
			if (t!=sf->hyperperiod) out<<",";
		}
		out<<"];"<<endl;
		out<<endl;
	}
}

void FileWriter::WriteRbfAndDbf2(Stateflow** sfs, int num, ostream& out) {
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		out<<"Stateflow-"<<i<<endl;

		// write rbf(t)
		out<<"%=================rbf(t)======================"<<endl;
		for (int t=0; t <= sf->hyperperiod; t+=sf->gcd) {
			out<<"rbf("<<t<<")="<<sf->get_rbf(t)<<endl;
		}
		out<<endl;

		// write dbf(t)
		out<<"%=================dbf(t)======================"<<endl;
		for (int t=0; t <= sf->hyperperiod; t+=sf->gcd) {
			out<<"dbf("<<t<<")="<<sf->get_dbf(t)<<endl;
		}
		out<<endl;
	}
}

void FileWriter::WriteResponseTimes(Stateflow** sfs, int num, ostream& out) {
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		out<<"Write the response times of Stateflow-"<<i<<endl;
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* transition = *iter;
			out << transition->src->name << " -> " << transition->snk->name 
				<< " [label=\" " << transition->wcet << " / " << transition->period << " \"]" << endl;
			out<<"Response times by rbf with static offsets"<<endl;
			for (map<int,double>::iterator m_iter = transition->responseTimeByRbfWithStaticOffsets.begin(); m_iter != transition->responseTimeByRbfWithStaticOffsets.end(); m_iter++)
				out<<m_iter->first<<" => "<<m_iter->second<<endl;
			out<<"Response times by rbf with arbitrary offsets"<<endl;
			out<<transition->responseTimeByRbfWithArbitraryOffsets<<endl;

			out<<"Response times by ibf with static offsets"<<endl;
			for (map<int,double>::iterator m_iter = transition->responseTimeByIbfWithStaticOffsets.begin(); m_iter != transition->responseTimeByIbfWithStaticOffsets.end(); m_iter++)
				out<<m_iter->first<<" => "<<m_iter->second<<endl;
			out<<"Response times by ibf with arbitrary offsets"<<endl;
			out<<transition->responseTimeByIbfWithArbitraryOffsets<<endl;
			out<<endl;
		}
		out<<endl;
	}
	out<<endl;
}