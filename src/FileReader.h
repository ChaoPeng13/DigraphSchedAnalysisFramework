/* \file FileReader.h
 * Generate digraphs or stateflows from files
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 11-sept-2015 : initial revision (CP)
 *
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_

//#include <vld.h>

#include <iostream>
#include <fstream>

#include "Digraph.h"
#include "Stateflow.h"

class FileReader {
public:
	static void DotFileReader(Digraph** &digraphs, int &num, const char* fname);
	static void DotFileReader(Stateflow** &sfs, int &num, int scale, const char* fname);
};

#endif