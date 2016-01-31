/* \file FileWriter.h
 * Save digraphs or stateflows into files
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 11-sept-2015 : initial revision (CP)
 *
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include <iostream>
#include <fstream>

#include "Digraph.h"
#include "Stateflow.h"

class FileWriter {
public:
	static void DotFileWriter(Digraph** digraphs, int num, const char* fname);
	static void DotFileWriter(Stateflow* sf, const char* fname);
	static void DotFileWriter(Stateflow** sfs, int num, const char* fname);
	/// \biref write a .m file to draw line picture
	static void DotFileWriter(Stateflow* sf, int range, string prefix);

	/// \brief write rbf and dbf
	static void WriteRbfAndDbf(Stateflow** sfs, int num, ostream& out);
	static void WriteRbfAndDbf2(Stateflow** sfs, int num, ostream& out);

	/// \brief write the response times
	static void WriteResponseTimes(Stateflow** sfs, int num, ostream& out);
};

#endif