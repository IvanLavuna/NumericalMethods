//
// Created by ivan on 16.10.20.
//

#ifndef NUMERICALMETHODSLABS_BASIC_H
#define NUMERICALMETHODSLABS_BASIC_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <initializer_list>
// TODO: after writing Matrix class integrate this with this code
//#include "Matrix.h"


/// Some defines to work more easily with flags
#define SET_UNDEFINED_BEHAVIOUR _funcState = UNDEFINED_BEHAVIOUR
#define SET_ONE_SOLUTION _funcState		  = ONE_SOLUTION
#define SET_INFINITY_SOLUTIONS _funcState  = INFINITY_SOLUTIONS
#define SET_NO_SOLUTION _funcState 		  = NO_SOLUTION

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define ERROR_MESSAGE(DESC) ("error: " AT "\n" "Description: " DESC "\n")



namespace nm
{
	using std::cin;
	using std::cout;
	using std::endl;
	using std::pair;
	using std::abs;
	using std::vector;

/** simple abstractions for better understanding code**/
	typedef vector<vector<double>> Matrix;
	typedef vector<double> RowMatrix;
	typedef vector<std::complex<double>> ComplexRow;
	typedef vector<vector<std::complex<double>>> ComplexMatrix;
	typedef long long ll;
	/** structure which identifies P, L and U matrices in
 * PA = LU factorization **/

	struct P_L_U
	{
		Matrix P;
		Matrix L;
		Matrix U;
	};
	struct L_U
	{
		Matrix L;
		Matrix U;
	};

/** some constant variables **/
	const double EPS = 1e-9;
	const int INF = 10000000;

/** Flags to work with functions **/
/** Depending on what type of function you call
 * different flag will be set or unset
 * **/

	enum FLAG
	{
		ONE_SOLUTION = 0,
		INFINITY_SOLUTIONS,
		UNDEFINED_BEHAVIOUR,
		NO_SOLUTION
	}; /// variable that defines state of function

	/// global variable
	extern FLAG _funcState;


/** helpful functions **/
	void printMatrix(Matrix &m);

	void printRowMatrix(RowMatrix &m);

	void printComplexRowMatrix(ComplexRow &m);
/**
 * Returns true if matrix is summetric otherwise returns false
 * **/
	bool isSymmetric(Matrix &A);

/**
 * Transparent matrix A is matrix whose rows equals to columns of
 * same matrix A(swapping rows with columns)
 * **/
	Matrix getTransparentMatrix(Matrix &A);

	Matrix getTransparentMatrix(RowMatrix &A);

/**
 *
 * @param A -> matrix that will be copied
 * @param n -> copy first n rows
 * @param m -> copy first m columns
 * @return matrix A[0 : n - 1][0 : m - 1]
 */
	Matrix copyPart(Matrix &A, int n, int m);

}

#endif //NUMERICALMETHODSLABS_BASIC_H
