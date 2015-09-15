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
		cout<<endl;
	}
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