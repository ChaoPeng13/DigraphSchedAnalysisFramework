/* \file DigraphExample.h
 * this file implements some digraph examples
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 16-sept-2015 : initial revision (CP)
 *
 */

#ifndef DIGRAPHEXAMPLE_H_
#define DIGRAPHEXAMPLE_H_

#include "Digraph.h"

class DigraphExample {
public:
	/**
	 * Used to test deep first search
	 * The digraph comes from Figure 22-4 of "Cormen et al: Introduction to
	 * agorithms", Chapter 22.3.
	 */
	static Digraph* generateDigraph0();
	/**
	 * Used to test strongly connected components algortihm
	 * The digraph comes from Figure 22-9 of "Cormen et al: Introduction to
	 * agorithms", Chapter 22.5. 
	 */
	static Digraph* generateDigraph1();
	/**
	 * An example of the paper RTNS2015.
	 * Used to calculate the linear factor (or utilization)
	 */
	static Digraph* generateDigraph2();
	/**
	 * An example of the paper RTS2014.
	 * Used to test rbf and dbf
	 */
	static Digraph* generateDigraph3();
	/**
	 * An example of Martin Stigge et al., The digraph real-time task model, RTAS2011
	 * Used to dbf
	 */
	static Digraph* generateDigraph4();
	/**
	 * Figure 5 in the paper "Schedulability analysis of periodic tasks implementing synchronous finite state machines"
	 * Used to test rbf
	 */
	static Digraph* generateDigraph5();
	/**
	 * A strongly connected example for 
	 * Figure 5 in the paper "Schedulability analysis of periodic tasks implementing synchronous finite state machines"
	 * by adding some loose edges
	 */
	static Digraph* generateDigraph6();
	 
};

#endif