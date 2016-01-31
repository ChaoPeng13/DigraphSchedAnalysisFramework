#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "MaxPlusAlgebra.h"
#include "Utility.h"

// define the positive infinity of double
extern double POS_INFINITY; // = std::numeric_limits<double>::infinity();
// define the negative infinity of double
extern double NEG_INFINITY; // = - std::numeric_limits<double>::infinity();

//extern bool compare_tow_matrices(double** A, double** B, int n);

/**
 * Test an implementation of Theorem 3.6 in M. Gavalec, linear matrix period in max-plus algebra,�Linear Algebra and its Applications, 
 * 307(1-3):167�82, 2000.
 * Example 2 in the paper.
 */
TEST(MaxPlusAlgebraTest, LinearMatrixPeriod)
{
	double** A = new double*[6];
	for (int i=0; i<6; i++) A[i] = new double[6];
	for (int i=0; i<6; i++) for (int j=0; j<6; j++) A[i][j] = NEG_INFINITY;

	A[0][1] = 3;
	A[1][2] = 0;
	A[1][4] = 2;
	A[2][0] = -1;
	A[2][3] = -3;
	A[3][0] = 4;
	A[4][5] = -4;
	A[5][3] = -1;
	A[5][4] = 6;

	int lfac = 1;
	int lper = MaxPlusAlgebra::calculate_linear_period(A,6,1,false);
	EXPECT_EQ(lper,4);
	
	std::map<int,double**> matrices;
	matrices[1]=A;
	int ldef = MaxPlusAlgebra::calculate_linear_defect(matrices,6,1,lper,100);
	//std::cout<<"ldef="<<ldef<<std::endl;
	
	for (int t=2; t<=100; t++) {
		matrices[t] = MaxPlusAlgebra::multiply_maxplus_matrix(matrices[(t-1)],matrices[1],6);
	}

	for (int t=ldef+1; t<=50; t++) {
		//std::cout<<"=================="<<t<<"====================="<<std::endl;

		double** temp = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(matrices[t],6,lfac,lper);
		double** temp2 = matrices[t+lper]; 
		
		EXPECT_EQ(Utility::compare_two_matrices(temp,temp2,6), true);
	}
	
}

/**
 * Test linear matrix defect
 * An example comes from Zeng and Di Natale, Computing periodic request functions to speed-up the analysis of non-cyclic task models.
 * RTS2014
 */
TEST(MaxPlusAlgebraTest, LinearMatrixDefect)
{
	double** A = new double*[5];
	for (int i=0; i<5; i++) A[i] = new double[5];
	for (int i=0; i<5; i++) for (int j=0; j<5; j++) A[i][j] = NEG_INFINITY;

	A[0][0] = 0.1;
	A[0][3] = 0.1;
	A[1][1] = 0;
	A[1][2] = 0.2;
	A[2][0] = 0.1;
	A[2][2] = 0;
	A[2][4] = 0.1;
	A[3][1] = 0;
	A[4][1] = 0;

	double lfac = 0.1;
	int lper = MaxPlusAlgebra::calculate_linear_period(A,5,lfac,false);
	//std::cout<<"lper="<<lper<<std::endl;
	EXPECT_EQ(lper,1);

	std::map<int,double**> matrices;
	matrices[1]=A;
	int ldef = MaxPlusAlgebra::calculate_linear_defect(matrices,5,lfac,lper,100);
	//Utility::output_matrix(matrices[6],5,5);
	EXPECT_EQ(ldef,6);
	
	for (int t=ldef; t<=10; t++) {
		//std::cout<<"=================="<<t-1<<"====================="<<std::endl;

		double** temp = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(matrices[t],5,lfac,lper);
		double** temp2 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices,t+lper,5); 
		double** temp3 = MaxPlusAlgebra::multiply_maxplus_matrix(matrices[1],matrices[t],5);
		//Utility::output_matrix(temp,5,5);
		//Utility::output_matrix(temp2,5,5);
		//Utility::output_matrix(temp3,5,5);
		//Utility::output_matrix(matrices[t],5,5);
		
		EXPECT_EQ(Utility::compare_two_matrices(temp,temp2,5), true);
	}
	
}

#endif