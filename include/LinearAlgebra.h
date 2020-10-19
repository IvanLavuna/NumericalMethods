//
// Created by Ivan on 06.10.20.
//
/** NOTE! All functions are not normally tested so bugs may occur **/

// TODO : Test normally each function


#ifndef NUMERICALMETHODSLABS_LAB1_H
#define NUMERICALMETHODSLABS_LAB1_H

#include "Basic.h"

/** Numerical methods space **/
namespace nm
{


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
	// TODO : insert this function to Basic.h
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

/**
 * Підраховує число обумовленості матриці
 * Чим більшим є число cond A , ти більша може бути відносна похибка
 * cond(A) = ||A|| * inverse(||A||)
 **/
 // TODO: insert this function into Basic.h
	double conditionNumber(Matrix& A);

	/**
	 * Matrix norm is the biggest sum in all rows
	 * **/
	// TODO: insert this function into Basic.h
	double MatrixNorm(Matrix& A);

	/**
	 * vector norm of order p is
	 * sqrt( sum( (x[i]) ^ p) ^ (1/p)), i = 1 ... n
	 * maks norm = max(abs(x[i])), i = 1 ... n
	 * @param A
	 * @return
	 */
	// TODO: insert this function into Basic.h
	double VectorNorm(RowMatrix& A,int p);

	/**
	 * Solves SLAE with complex numbers using Gauss-Jordan
	 **/
	ComplexRow solveSLAEWithComplexElements(ComplexMatrix A);
}



#endif //NUMERICALMETHODSLABS_LAB1_H











