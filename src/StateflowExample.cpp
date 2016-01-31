#include "StateflowExample.h"

/**
 * An example for Figure 1 in 
 * Zeng and Di Natale, Schedulability analysis of periodic tasks implementing synchronous finite state machines.
 * ECRTS2012.
 */
Stateflow* StateflowExample::generateStateflow0() {
	Stateflow* sf = new Stateflow(100);

	State* s1 = new State("s1",1);
	State* s2 = new State("s2",2);
	State* s3 = new State("s3",3);

	sf->add_state(s1);
	sf->add_state(s2);
	sf->add_state(s3);

	Transition* t12 = new Transition(s1,s2,"a1");
	t12->wcet = 25;
	t12->period = 200;
	t12->priority = 1;
	sf->add_transition(t12);

	Transition* t23_0 = new Transition(s2,s3,"a3");
	t23_0->wcet = 10;
	t23_0->period = 200;
	t23_0->priority = 1;
	sf->add_transition(t23_0);

	Transition* t23_1 = new Transition(s2,s3,"a4");
	t23_1->wcet = 15;
	t23_1->period = 500;
	t23_1->priority = 2;
	sf->add_transition(t23_1);

	Transition* t31 = new Transition(s3,s1,"a2");
	t31->wcet = 30;
	t31->period = 500;
	t31->priority = 2;
	sf->add_transition(t31);

	if (false) sf->write_graphviz(cout);
	return sf;
}

Stateflow* StateflowExample::generateStateflow1() {
	Stateflow* sf = new Stateflow(100);

	State* s1 = new State("s1",1);
	sf->add_state(s1);

	Transition* t11 = new Transition(s1,s1);
	t11->wcet = 25;
	t11->period = 400;
	t11->priority = 1;
	sf->add_transition(t11);

	if (false) sf->write_graphviz(cout);

	return sf;
}

Stateflow* StateflowExample::generateStateflow2() {
	Stateflow* sf = new Stateflow(100);

	State* s1 = new State("s1",1);
	sf->add_state(s1);

	Transition* t11 = new Transition(s1,s1);
	t11->wcet = 150;
	t11->period = 200;
	t11->priority = 1;
	sf->add_transition(t11);

	return sf;
}

/**
 * A modeifed example for exchanging the wcets of a1 and a3 for Figure 1 in 
 * Zeng and Di Natale, Schedulability analysis of periodic tasks implementing synchronous finite state machines.
 * ECRTS2012.
 */
Stateflow* StateflowExample::generateStateflow3() {
	Stateflow* sf = new Stateflow(100);

	State* s1 = new State("s1",1);
	State* s2 = new State("s2",2);
	State* s3 = new State("s3",3);

	sf->add_state(s1);
	sf->add_state(s2);
	sf->add_state(s3);

	Transition* t12 = new Transition(s1,s2,"a1");
	t12->wcet = 10;
	t12->period = 200;
	t12->priority = 1;
	sf->add_transition(t12);

	Transition* t23_0 = new Transition(s2,s3,"a3");
	t23_0->wcet = 25;
	t23_0->period = 200;
	t23_0->priority = 1;
	sf->add_transition(t23_0);

	Transition* t23_1 = new Transition(s2,s3,"a4");
	t23_1->wcet = 15;
	t23_1->period = 500;
	t23_1->priority = 2;
	sf->add_transition(t23_1);

	Transition* t31 = new Transition(s3,s1,"a2");
	t31->wcet = 30;
	t31->period = 500;
	t31->priority = 2;
	sf->add_transition(t31);

	if (false) sf->write_graphviz(cout);
	return sf;
}