//
// Created by Ivan on 06.10.20.
//
/** NOTE! All functions are not normally tested so bugs may occur **/

// TODO : Test normally each function


#ifndef NUMERICALMETHODSLABS_LAB1_H
#define NUMERICALMETHODSLABS_LAB1_H




#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>


/// Some defines to work more easily with flags
#define SET_UNDEFINED_BEHAVIOUR state = UNDEFINED_BEHAVIOUR
#define SET_ONE_SOLUTION state 		  = ONE_SOLUTION
#define SET_INFINITY_SOLUTIONS state  = INFINITY_SOLUTIONS
#define SET_NO_SOLUTION state 		  = NO_SOLUTION

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define ERROR_MESSAGE(DESC) ("error: " AT "\n" "Description: " DESC "\n")


/** Numerical methods space **/
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

// TODO : rewrite every function with flags
	enum FLAG
	{
		ONE_SOLUTION = 0,
		INFINITY_SOLUTIONS,
		UNDEFINED_BEHAVIOUR,
		NO_SOLUTION
	} state; /// variable that defines state of function


/** helpful functions **/
	void printMatrix(Matrix &m);

	void printRowMatrix(RowMatrix &m);

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

///////////////////////////////////////////////////////////

/** LINEAR ALGEBRA **/
/**
 * @brief
 * 	solves SLAE of equation using Gauss method
 * @return
 * first - number of solutions
 * second - solution vector if number of solutions is one
 *
 * @param
 * matrix which defines SLAE
 */
	RowMatrix solveSLAEByGauss(Matrix matrix);


/**
 * Calculates determinant of quadratic matrix using Gauss method
 *	throws an exception of type invalid argument if matrix is not quadratic
 *
 **/
	double getDeterminant(Matrix matrix);


/**
 *
 * The rank of a matrix is the largest number of linearly independent
 * rows/columns of the matrix.
 * @tparam T
 * @param matrix
 * @return rank
 * Method also uses Gauss elimination method
 */
	int getRank(Matrix matrix);

/**
 * @brief
 * This function returns L matrix and U matrix.
 * This function is not using rows exchanging!
 * If module of some diagonal element is less that EPS,
 * it will return pair of zero matrices
 * Obviously, it works only on squared matrices
 */

	L_U LU_factorization(Matrix matrix);


/**
 * @brief
 * 		This function performs PA = LU factorization
 * 		(with swapping some rows)
 * 		and returns L and U matrices
 * @param
 * 		vector of size n with vectors of size n
 * @exception
 * 		throws an exception if matrix is not squared
 *
 * @return
 * 		returns L and U matrices if they exist, else
 * 		returns pair of zero matrices
 */
	L_U PA_LU_factorization(Matrix A);

/** returns P L and U matrices in PA = LU factorization **/
	P_L_U get_P_L_U(Matrix A);

/**
 * @brief
 * 		multiplies two matrices
 * @exception
 * 		returns exception if number of columns of
 * 		first matrix A is not equal to number of rows
 * 		in matrix B
 * 		time complexity is O(n^3)
 */
	Matrix multiply(Matrix &A, Matrix &B);

/**
 * @brief
 * 		Solves SLAE using PA = LU factorisation
 * 		Works only on quadratic matrices
 *
 * 	algorithm description
 * 		LUX = PB
 * 		LC  = PB
 * 		UX  = C
 * @exception
 * 		throws an exception if matrix A is not quadratic.
 * @param
 * 		Parameter is matrix which identifies A and B
 * 		matrices;
 * @returns
 * 	RowMatrix size n if solution exist and is single,
 * 	else return ZeroRowMatrix
 */
	RowMatrix solveSLAEByLUFactorization(Matrix &matrix);

/**
 *	finds inverse of a matrix using Gauss-Jordan
 *	algorithm.
 *	@returns zero matrix if there is not exist inverse matrix
 *	Gives @exception if matrix is not quadratic
 **/
	Matrix getInverseMatrix(Matrix A);


/**
 * @brief
 * 	Returns LU( A = L * L_tr) factorization by Cholesky
 *	Matrix must be symmetric. If factorization is not
 *	possible if returns  zero matrices
 *
 *	@Note
 *		Maybe for this algorithm to work elements of matrix must be positive
 *	Throws an @exception if matrix is not symmetric
 *
 *	@returns zero matrix if there is NOT exist LU factorization
 **/
	L_U get_L_U_factorizationByCholesky(Matrix A);

/**
 * @param is matrix M which defines SLAE
 *
 * @brief
 *	Solves SYMMETRIC SLAE by Cholesky
 *
 * @exception
 *	Throws an exception if matrix A is not symmetric
 *
 */
	RowMatrix solveSLAEByCholesky(Matrix &M);


/**
 * Розв'язання матричного рівняння за допомогою методу оборотів
 *
 */
	RowMatrix solveSLAEByMethodOfTurns(Matrix A);
}

#endif //NUMERICALMETHODSLABS_LAB1_H











