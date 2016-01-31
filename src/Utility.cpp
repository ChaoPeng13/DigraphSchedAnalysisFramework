#include "Utility.h"

double Utility::EPSILON = 0.000001;

// return the greatest common divisor for two integers
int Utility::math_gcd(int a, int b) {
		if (b==0) return a;
		return math_gcd(b, a%b);
}

int Utility::math_lcm(int a, int b) {
		if (b==0) return a;
		return (a*(b/math_gcd(a,b)));
}

std::string Utility::int_to_string(int a) {
		char temp[10];
		sprintf(temp,"%d",a);
		return temp;
}

double** Utility::creat_matrix(int nrow, int ncol) {
	double** A = new double*[nrow];
	for (int i=0; i<nrow; i++) A[i] = new double[ncol];
	return A;
}

void Utility::output_matrix(double** A, int nrow, int ncol) {
	for (int i=0; i<nrow; i++) {
		for (int j=0; j<ncol; j++) {
			std::cout<< A[i][j];
			if (j==ncol-1) std::cout<<std::endl;
			else std::cout<<"\t";
		}
	}
}

double* Utility::uniformly_distributed(int n, double tUtil) {
	double* ret = new double[n];

	double sum = tUtil;
	for (int i=0; i<n-1; i++) {
		double nextsum = sum*pow((double)rand()/RAND_MAX, 1.0/(n-i-1));
		ret[i] = sum-nextsum;
		sum = nextsum;
	}

	ret[n-1] = sum;
	return ret;
}

bool Utility::compare_two_matrices(double** A, double** B, int n) {
	for (int i=0; i<n; i++) for (int j=0; j<n; j++)
		if (abs(A[i][j]-B[i][j])>EPSILON) {
			if (false) {
				std::cout<<i<<","<<j<<std::endl;
				std::cout<<A[i][j]<<"..."<<B[i][j]<<std::endl;
			}
			return false;
		}
	return true;
}