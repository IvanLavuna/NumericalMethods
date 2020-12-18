//
// Created by ivan on 16.10.20.
//

#ifndef NUMERICALMETHODSLABS_BASIC_H
#define NUMERICALMETHODSLABS_BASIC_H

#include <complex>
#include "Matrix.h"

namespace nm
{
	using std::cin;
	using std::cout;
	using std::endl;
	using std::pair;
	using std::abs;
	using std::vector;

/** simple abstractions for better understanding code**/


	typedef vector<double> RowMatrix;
	typedef vector<vector<std::complex<double>>> ComplexMatrix;
	typedef vector<std::complex<double>> ComplexRow;
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
	bool isSymmetric(const Matrix &A);

/**
 * Transparent matrix A is matrix whose rows equals to columns of
 * same matrix A(swapping rows with columns)
 * **/
	Matrix TransparentMatrix(Matrix &A);

	Matrix TransparentMatrix(RowMatrix &A);

/**
 *
 * @param A -> matrix that will be copied
 * @param n -> copy first n rows
 * @param m -> copy first m columns
 * @return matrix A[0 : n - 1][0 : m - 1]
 */
	Matrix copyPart(Matrix &A, int n, int m);

	/**
	 * @param n - size of squared matrix
	 * @return Identity matrix of size n
	 */
	Matrix IdentityMatrix(int n);

}

#endif //NUMERICALMETHODSLABS_BASIC_H
