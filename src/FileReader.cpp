#include "FileReader.h"

void FileReader::DotFileReader(Digraph** &digraphs, int &num, int scale, const char* fname) {
	ifstream fin(fname);
	if (!fin.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}
	// some parameters
	string start("digraph G {");
	string end("}");
	string arrow("->");

	num = -1; // count of Digraphs
	int index = 0; // index of node
	vector<Digraph*> dg_vec;
	Digraph* digraph;
	bool newDG = false;

	map<string, Node*> name_to_node;

	string line;
	while (getline(fin,line)) { // read one line into string line
		if (!line.compare(start)) { // start one digraph
			newDG = true;
			num++;
			continue;
		}

		if (!line.compare(end)) // end one stateflow
			continue;

		if (line.length() == 0) { // empty line
			name_to_node.clear();
			index = 0;
			continue;
		}

		if (newDG) {
			digraph = new Digraph(num);
			dg_vec.push_back(digraph);
			newDG = false;
		}

		// parse the line
		stringstream ss(line);
		string sub_str;
		vector<string> str_list;
		while(getline(ss,sub_str,' ')) {
			str_list.push_back(sub_str);
		}

		int position = line.find(arrow);
		if (position != string::npos) { // read one edge
			Node* src = name_to_node[str_list.at(0)];
			Node* snk = name_to_node[str_list.at(2)];
			Edge* edge = new Edge(src, snk);
			edge->separationTime = atof(str_list.at(4).c_str());

			digraph->add_edge(edge);
		} else { // read one node
			string name = str_list.at(0);
			int wcet = atof(str_list.at(4).c_str());
			int deadline = atof(str_list.at(6).c_str());

			Node* node = new Node(name, index++, scale, wcet, deadline);
			name_to_node[name] = node;

			digraph->add_node(node);
		}
	}
	fin.close();
	num = dg_vec.size();
	digraphs = new Digraph*[num];
	for (int i=0; i<num; i++)
		digraphs[i] = dg_vec.at(i);
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

Stateflow* FileReader::ReadOneStateflow(int scale, const char* fname) {
	ifstream fin(fname);
	if (!fin.is_open()) {
		cerr << "Can't open "<<fname<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	// some parameters
	string start("digraph G {");
	string end("}");
	string arrow("->");

	int index = 0; // index of state
	Stateflow* sf = new Stateflow(scale);

	map<string, State*> name_to_state;

	string line;
	while(getline(fin,line)) { // read one line into string line
		if (!line.compare(start)) continue;
		if (!line.compare(end)) continue; 
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
	
	return sf;
}