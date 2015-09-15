#include "FileReader.h"

void FileReader::DotFileReader(Digraph** &digraphs, int &num, const char* fname) {
	ifstream fin(fname);
	if (!fin.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	
}

void FileReader::DotFileReader(Stateflow** &sfs, int &num, int scale, const char* fname) {
	ifstream fin(fname);
	if (!fin.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	// some parameters
	string start("digraph G {");
	string end("}");
	string arrow("->");

	num = -1; // count of Stateflows
	int index = 0; // index of state
	vector<Stateflow*> sf_vec;
	Stateflow* sf;
	bool newSF = false;

	map<string, State*> name_to_state;

	string line;
	while(getline(fin,line)) { // read one line into string line
		if (!line.compare(start)) { // start one stateflow
			newSF = true;
			num++;
			continue;
		}

		if (!line.compare(end))  // end one stateflow
			continue;

		if (line.length() == 0) { // empty line
			name_to_state.clear();
			index = 0;
			continue;
		}

		if (newSF) {
			sf = new Stateflow(num,scale);
			sf_vec.push_back(sf);
			newSF = false;
		}

		// parse the line
		stringstream ss(line);
		string sub_str;
		vector<string> str_list;
		while(getline(ss,sub_str,' ')) {
			str_list.push_back(sub_str);
		}

		int position = line.find(arrow);
		if (position != string::npos) { // Read one transition
			State* src = name_to_state[str_list.at(0)];
			State* snk = name_to_state[str_list.at(2)];
			Transition* t = new Transition(src, snk);
			t->wcet = atof(str_list.at(4).c_str());
			t->period = atof(str_list.at(6).c_str())*scale;

			sf->add_transition(t);
		} else { // Read one state
			string name = str_list.at(0);
			State* state = new State(name, index++);
			name_to_state[name] = state;

			sf->add_state(state);
		}
	}

	fin.close();
	
	num = sf_vec.size();
	sfs = new Stateflow*[num];
	for (int i=0; i<num; i++)
		sfs[i] = sf_vec.at(i);
}