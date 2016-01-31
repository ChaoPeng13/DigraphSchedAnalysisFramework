/* \file StateflowExample.h
 * this file implements some Stateflow examples
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 16-sept-2015 : initial revision (CP)
 *
 */

#ifndef STATEFLOWEXAMPLE_H_
#define STATEFLOWEXAMPLE_H_

#include "Stateflow.h"

class StateflowExample {
public:
	/**
	 * An example for Figure 1 in 
	 * Zeng and Di Natale, Schedulability analysis of periodic tasks implementing synchronous finite state machines.
	 * ECRTS2012.
	 */
	 static Stateflow* generateStateflow0();

	 static Stateflow* generateStateflow1();
	 static Stateflow* generateStateflow2();

	 static Stateflow* generateStateflow3();

};

#endif