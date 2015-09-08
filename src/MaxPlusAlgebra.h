/* \file MaxPlusAlgebra.h
*  this file implements some algorithms for linear matrix. 
*  \author Chao Peng
*  
*  Changes
*  ------
*  25-aug-2015 : initial revision (CP)
*
*/

#ifndef MAXPLUSALGEBRA_H
#define MAXPLUSALGEBRA_H

#include <algorithm>
#include <vector>
#include <map>
#include <set>

class MaxPlusAlgebra {
public:
	static double** calculate_metric_matrix(double** B, int n);
	static int calculate_linear_period(double** A, int n, double lfac, bool debug);
	static int calculate_gcd_cycle_length_hcc(double** A, int n, std::set<int> hcc, bool debug);
	static int calculate_linear_defect(std::map<int, double**>& m_map, std::map<int, int> &me_map, int n, double lfac, int lper, int tf);
	static int calculate_linear_defect(std::map<int, double**>& m_map, int n, double lfac, int lper, int tf);
	static double maximum_element(double** A, int n);
	static double** periodicly_calculate_maxplus_matrix_power(double** A, int n, double lfac, int lper);
	static double** multiply_maxplus_matrix(double** A, double** B, int n);
	static double** multiply_maxplus_matrix(std::map<int,double**> & matrices, int t, int n);
};

#endif