#include "MaxPlusAlgebra.h"

#include "Utility.h"
#include <queue>

// define the positive infinity of double
double POS_INFINITY = std::numeric_limits<double>::infinity();
// define the negative infinity of double
double NEG_INFINITY = - std::numeric_limits<double>::infinity();
// define a very small double
double EPSILON = 0.000001;


double** MaxPlusAlgebra::calculate_metric_matrix(double** B, int n) {
	double** barB = new double*[n];
	for (int i=0; i<n; i++) barB[i] = new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) barB[i][j] = B[i][j];
	for (int k=0; k<n; k++) for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
		//std::cout<< barB[i][j] <<"\t"<< barB[i][k] <<"\t"<<barB[k][j]<<std::endl;
		//double a = barB[i][j];
		//double b = barB[i][k];
		//double c = barB[k][j];
		barB[i][j] = std::max(barB[i][j], barB[i][k]+barB[k][j]);
		//a = barB[i][j];
		//b = barB[i][k]; 
		//std::cout<< barB[i][j] << std::endl;
	}

	return barB;
}

// calculating linear period:
// an implementation of Theorem 3.6 in M. Gavalec, linear matrix period in max-plus algebra,�Linear Algebra and its Applications, 
// 307(1-3):167�82, 2000.
// We should operate on the execution request matrix whithout zero-cycle
int MaxPlusAlgebra::calculate_linear_period (double** A, int n, double lfac, bool debug) {
	if (debug) {
		std::cout<< "A="<<std::endl;
		Utility::output_matrix(A,n,n);
	}
	// step 1: B=A-\lambda(A)
	double** B = new double*[n];
	for (int i=0; i<n; i++) B[i]=new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) B[i][j] = A[i][j]-lfac;

	if(debug) {
		std::cout<< "B=" <<std::endl;
		Utility::output_matrix(B,n,n);
	}

	// step 2: metric matrix of B, barB[i][j] is the maximal weight of a walk from i to j
	double** barB = calculate_metric_matrix(B,n);
	if(debug) {
		std::cout<< "barB=" <<std::endl;
		Utility::output_matrix(barB,n,n);
	}

	// step 3: the highly connected components of the digraph coresponding to B
	double** C = new double*[n];
	for (int i=0; i<n; i++) C[i]=new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) C[i][j]= barB[i][j]+barB[j][i];

	if(debug) {
		std::cout<< "C=" <<std::endl;
		Utility::output_matrix(C,n,n);
	}

	bool* isVisited = new bool[n];
	for (int i=0; i<n; i++) isVisited[i] = false;
	std::vector<std::set<int>> hcc_vec;
	for (int i=0; i<n; i++) {
		bool trivial = true;
		if (isVisited[i]) continue;
		std::set<int> hcc;

		for (int j=0; j<n; j++) {
			if ((int)(C[i][j]*1000) == 0) {
				trivial = false;
				break;
			}
		}

		if (!isVisited[i] && !trivial) {
			hcc.insert(i);
			isVisited[i] = true;
			for (int j=0; j<n; j++) {
				if ((int)(C[i][j]*1000) == 0) {
					hcc.insert(j);
					isVisited[j] = true;
				}
			}
		}
		hcc_vec.push_back(hcc);
	}

	// setp 4: find no zero-cycle in G(B)
	double** D = new double*[n];
	for (int i=0; i<n; i++) D[i]=new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) D[i][j] = B[i][j]+barB[j][i];

	if(debug) {
		std::cout<< "D=" <<std::endl;
		Utility::output_matrix(D,n,n);
	}
	
	// step 5: linear period, lcm{hper(hcc)}
	double** dotB = new double*[n];
	for (int i=0; i<n; i++) dotB[i]=new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
		if (D[i][j]<0 && D[i][j]!=NEG_INFINITY) dotB[i][j] = NEG_INFINITY;
		else dotB[i][j] = B[i][j];
	}

	if(debug) {
		std::cout<< "dotB=" <<std::endl;
		Utility::output_matrix(dotB,n,n);
	}

	// find all cycles with maximum cycle mean for highly connected components
	int lcm = 1;
	for (std::vector<std::set<int>>::iterator hcc_iter = hcc_vec.begin(); hcc_iter != hcc_vec.end(); hcc_iter++) {
		std::set<int> hcc = *hcc_iter;
		int gcd = calculate_gcd_cycle_length_hcc(dotB,n,hcc,debug);
		//std::cout<<"gcd = "<<gcd<<std::endl;
		lcm = Utility::math_lcm(lcm, gcd);
	}

	// release memory
	// release B
	for (int i=0; i<n; i++) delete[] B[i];
	delete[] B;
	// release barB
	for (int i=0; i<n; i++) delete[] barB[i];
	delete[] barB;
	// release C
	for (int i=0; i<n; i++) delete[] C[i];
	delete[] C;
	// release isVisted
	delete[] isVisited;
	// release hcc_vec
	hcc_vec.clear();
	// release D
	for (int i=0; i<n; i++) delete[] D[i];
	delete[] D;
	// release dotB
	for (int i=0; i<n; i++) delete[] dotB[i];
	delete[] dotB;

	return lcm;
}

// Calculate the greatest common divisor of cycle lengths with \lambda(A) for highly connected components
int MaxPlusAlgebra::calculate_gcd_cycle_length_hcc(double** A, int n, std::set<int> hcc, bool debug) {
	std::vector<std::vector<int>> cycles;
	std::vector<std::vector<std::vector<int>>> omegas;

	// paths with 1 length
	std::vector<std::vector<int>> omega0;
	for (std::set<int>::iterator iter=hcc.begin(); iter != hcc.end(); iter++) {
		std::vector<int> path;
		path.push_back(*iter);
		omega0.push_back(path);
	}
	omegas.push_back(omega0);
	
	for(int i = 0; i<hcc.size(); i++) {
		std::vector<std::vector<int>> omega = omegas.back();
		std::vector<std::vector<int>> new_omega;
		for (std::set<int>::iterator iter = hcc.begin(); iter != hcc.end(); iter++) {
			int next = *iter;
			for (std::vector<std::vector<int>>::iterator iter1 = omega.begin(); iter1 != omega.end(); iter1++) {
				std::vector<int> path = *iter1;
				int last = path.back();
				if (A[last][next] != NEG_INFINITY) {
					bool isCycle = false;
					std::vector<int> new_path;
					for (std::vector<int>::iterator iter2 = path.begin(); iter2 != path.end(); iter2++) {
						new_path.push_back(*iter2);
						if (iter2 != path.begin()) {
							if (next == *iter2) isCycle = true;
						}
					}
					new_path.push_back(next);
					
					if (path.front() == next)
						cycles.push_back(new_path);
					else {
						if (!isCycle) new_omega.push_back(new_path);
					}
				}
			}
		}
		omegas.push_back(new_omega);
	}

	int gcd = 0;
	for (std::vector<std::vector<int>>::iterator iter = cycles.begin(); iter != cycles.end(); iter++) {
		std::vector<int> cycle = *iter;
		if (debug) {
			std::cout<<"cycle length="<<cycle.size()<<std::endl;
			for (std::vector<int>::iterator iter1 = cycle.begin(); iter1 != cycle.end(); iter1++) {
				std::cout<<*iter1<<"->";
			}
			std::cout<<std::endl;
		}
		gcd = Utility::math_gcd(gcd, cycle.size()-1);
	}
	
	return gcd;
}


/// \brief calculate the linear defect for a reducible matrix
/// Find the first r that satisfies A^(r+p） = A^(r) + p \times q
int MaxPlusAlgebra::calculate_linear_defect(std::map<int, double**> &m_map,  std::map<int, int> &me_map, int n, double lfac, int lper, int tf, int gcd) {
	double** A;
	for (int t=2; t<=tf; t++) {
		if (m_map.find(t) != m_map.end())
			A = m_map[t];
		else
		    A = MaxPlusAlgebra::multiply_maxplus_matrix(m_map[t-1],m_map[1],n);
		m_map[t] = A;
		double maximum_element = MaxPlusAlgebra::maximum_element(A, n);
		int temp = static_cast<int> (maximum_element);
		me_map[t] = temp;

		if (t>lper) {
			bool found = true;
			//std::cout<<"=============="<<t<<"=============="<<std::endl;
			double** B = m_map[t-lper];
			for ( int i=0; i<n; i++) {
				for (int j=0; j<n; j++) {
					double a = A[i][j];
					double b = B[i][j]+lfac*gcd*lper;
				
					//std::cout<<"a="<<a<<"\tb="<<b<<std::endl;
					if (abs(a-b)>EPSILON) {
						found = false;
						break;
					}
				}
				if (!found) break;
			}
			if(found) return t-lper;
		}
	}
	return -1; // used to identify there is no linear defect between (0,tf]
}

/// \brief calculate the linear defect for a reducible matrix
/// Find the first r that satisfies A^(r+p） = A^(r)  + p \times q
int MaxPlusAlgebra::calculate_linear_defect(std::map<int, double**> &m_map, int n, double lfac, int lper, int tf) {
	//std::cout<<"Calculate linear defect"<<std::endl;
	double** A;
	
	for (int t=2; t<=tf; t++) {
		if (m_map.find(t) != m_map.end())
			A = m_map[t];
		else
		    A = MaxPlusAlgebra::multiply_maxplus_matrix(m_map[t-1],m_map[1],n);
		m_map[t] = A;

		if (t>lper) {
			bool found = true;
			//std::cout<<"=============="<<t<<"=============="<<std::endl;
			double** B = m_map[t-lper];
			for ( int i=0; i<n; i++) {
				for (int j=0; j<n; j++) {
					double a = A[i][j];
					double b = B[i][j]+lfac*lper;
				
					//std::cout<<"a="<<a<<"\tb="<<b<<std::endl;
					if (abs(a-b)>EPSILON) {
						found = false;
						break;
					}
				}
			}
			if(found) return t-lper;
		}
	}

	return -1; // used to identify there is no linear defect between (0,tf]
}

double MaxPlusAlgebra::maximum_element(double** A, int n) {
	double max_element = NEG_INFINITY;
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) max_element = std::max(max_element, A[i][j]);

	return max_element;
}

double** MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(double** A, int n, double lfac, int lper) {
	double** B = new double*[n];
	for (int i=0; i<n; i++) B[i] = new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) B[i][j] = A[i][j] + lfac*lper;

	return B;
}

double** MaxPlusAlgebra::multiply_maxplus_matrix(double** A, double** B, int n) {
	double** C = new double*[n];
	for (int i=0; i<n; i++) C[i] = new double[n];

	for (int i=0; i<n; i++) for (int j=0; j<n; j++) C[i][j] = NEG_INFINITY;

	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++) 
			for (int k=0; k<n; k++)
				C[i][j] = std::max(C[i][j],A[i][k]+B[k][j]);

	return C;
}

double** MaxPlusAlgebra::power_maxplus_matrix(double** A, int n, int power) {
	double** C = new double*[n];
	for (int i=0; i<n; i++) C[i] = new double[n];

	for (int i=0; i<n; i++) for (int j=0; j<n; j++) C[i][j] = NEG_INFINITY;

	double** B = new double*[n];
	for (int i=0; i<n; i++) B[i] = new double[n];
	for (int i=0; i<n; i++) for (int j=0; j<n; j++) B[i][j] = A[i][j];

	for (int tp = 2; tp <= power; tp++) {
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++) 
				for (int k=0; k<n; k++)
					C[i][j] = std::max(C[i][j],A[i][k]+B[k][j]);

		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++) {
				B[i][j] = C[i][j];
				C[i][j] = NEG_INFINITY;
			}
	}

	for (int i=0; i<n; i++) delete[] B[i];
	delete[] B;

	return C;
}

double** MaxPlusAlgebra::multiply_maxplus_matrix(std::map<int,double**> & matrices, int t, int n) {
	if (matrices.find(t) != matrices.end())
		return matrices[t];

	double** C = new double*[n];
	for (int i=0; i<n; i++) C[i] = new double[n];

	for (int i=0; i<n; i++) for (int j=0; j<n; j++) C[i][j] = NEG_INFINITY;

	double** A = multiply_maxplus_matrix(matrices, t-1, n);
	double** B = matrices[1];

	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++) 
			for (int k=0; k<n; k++)
				C[i][j] = std::max(C[i][j],A[i][k]+B[k][j]);
	matrices[t] = C;
	return C;
}
