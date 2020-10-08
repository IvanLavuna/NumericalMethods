//
// Created by ivan on 06.10.20.
//

#ifndef NUMERICALMETHODSLABS_LAB1_H
#define NUMERICALMETHODSLABS_LAB1_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::pair;
using std::abs;
using std::vector;
typedef vector<vector<double>> Matrix;
const double EPS = 1e-9;
const int INF = 10000000;


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
Matrix Gauss(Matrix matrix)
{
	int n = matrix.size();
	int m = static_cast<int>(matrix[0].size()) - 1;

	vector<double> solution(m + 1);

	vector<bool> arb(m,false);

	for (int i = 0; i < n; ++i)
	{
		int pivot = i;
		/** searching for pivot element **/
		for (int j = i; j < n; ++j)
		{
			if (abs(matrix[pivot][i]) < abs(matrix[j][i]))
				pivot = j;
		}

		if(abs(matrix[pivot][i]) < EPS)
			continue;

		arb[i] = true;

		matrix[i].swap(matrix[pivot]);

		double div = matrix[i][i];

		for (int j = i; j <= m; ++j)
		{
			matrix[i][j] /= div;
		}

		for (int j = i + 1; j < n; ++j)
		{
			double mul = matrix[j][i];
			for (int k = 0; k <= m; ++k)
			{
				matrix[j][k] -= (mul * matrix[i][k]);
			}
		}
	}

	int ind = n - 1;

	solution[ind--] = matrix[n - 1][m];

	for (int i = n - 2; i >= 0; --i)
	{
		double cur_sol = matrix[i][m];
		for (int j = m - 1; j > i; --j)
		{
			cur_sol -= (matrix[i][j] * solution[j]);
		}
		solution[ind--] = cur_sol;
	}

	/** if mistake of some answer is bigger than some
	 * EPS , than there are not any solution**/
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		for (int j = 0; j < m; ++j)
			sum += solution[j] * matrix[i][j];
		if (abs (sum - matrix[i][m]) > EPS)
			return {0,solution};
	}

	/** checking for infinitely mane solutions **/
	for (int i = 0; i < m; ++i)
	{
		if(!arb[i]) return {INF,solution};
	}
	return {1,solution};
}


/**
 * Calculates determinant of quadratic matrix using Gauss method
 *	throws an exception of type invalid argument if matrix is not quadratic
 *
 **/
double Determinant(Matrix matrix)
{
	if(matrix.size() != matrix[0].size())
	{
		throw std::invalid_argument("Matrix must be quadratic\n");
	}

	int n = matrix.size();

	bool sign = true;
	double det = 1;

	for (int i = 0; i < n; ++i)
	{
		int pivot = i;
		/** searching for pivot element **/
		for (int j = i; j < n; ++j)
		{
			if (abs(matrix[pivot][i]) < abs(matrix[j][i]))
				pivot = j;
		}

		/** if pivot element is zero determinant is also zero **/
		if(std::abs(matrix[pivot][i]) < EPS)
			return 0;

		if(pivot != i)
		{
			sign = !sign;
			matrix[i].swap(matrix[pivot]);
		}

		/** multiplying determinant by some diagonal elem */
		det *= matrix[i][i];

		for (int j = i + 1; j < n; ++j)
		{
			double mul = matrix[j][i] / matrix[i][i];
			for (int k = 0; k < n; ++k)
			{
				matrix[j][k] -= (matrix[i][k] * mul);
			}
		}
	}
	return sign ? det : -det;
}


/**
 *
 * The rank of a matrix is the largest number of linearly independent
 * rows/columns of the matrix.
 * @tparam T
 * @param matrix
 * @return rank
 * Method also uses Gauss elimination method
 */


int Rank(Matrix matrix)
{
	int n = matrix.size();
	int m = static_cast<int>(matrix[0].size()) - 1;
	int rank = 0;

	for (int i = 0; i < n; ++i)
	{
		int pivot = i;
		for (int j = i; j < n; ++j)
		{
			if (abs(matrix[pivot][i]) < abs(matrix[j][i]))
				pivot = j;
		}

		if(std::abs(matrix[pivot][i]) < EPS)
			continue;

		/** incrementing runk **/
		rank++;

		/** swapping two rows**/
		matrix[i].swap(matrix[pivot]);

		double div = matrix[i][i];
		for (int j = i; j <= m; ++j)
		{
			matrix[i][j] /= div;
		}

		for (int j = i+1; j < n; ++j)
		{
			double mul = matrix[j][i];
			for (int k = 0; k <= m; ++k)
			{
				matrix[j][k] -= (mul * matrix[i][k]);
			}
		}
	}

	return rank;
}

/**
 * @brief
 * This function returns L matrix and U matrix.
 * This function is not using rows exchanging!
 * If module of some diagonal element is less that EPS,
 * it will return pair of zero matrices
 * Obviously, it works only on squared matrices
 */

pair<Matrix,Matrix> LU_factorization(Matrix matrix)
{
	if(matrix.size() != matrix[0].size())
	{
		throw std::invalid_argument("Matrix must be squared!");
	}

	int n = matrix.size();

	/// creating matrix already filled with zeros
	Matrix L(n,vector<double>(n,0) );

	/// filling diagonals with ones
	for (int i = 0; i < n; ++i) L[i][i] = 1;

	/// starting particularly like in Gauss elimination algorithm
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			if(abs(matrix[i][i]) < EPS)
			{
				return {Matrix (),Matrix ()};
			}
			double mul = matrix[j][i] / matrix[i][i];
			L[j][i] = mul;
			for (int k = 0; k < n; ++k)
			{
				matrix[j][k] -= (matrix[i][k] * mul);
			}
		}
	}
	return {L,matrix};
}


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
pair<Matrix ,Matrix> PA_LU_factorization(Matrix A)
{
	if(A.size() != A[0].size())
	{
		throw std::invalid_argument("Matrix must be squared!");
	}
	int n = A.size();

	/** creates and fills matrix with zeros**/
	Matrix L(n,vector<double>(n,0));


	/// starting particularly like in Gauss elimination algorithm
	for (int i = 0; i < n; ++i)
	{
		int pivot = i;
		for (int j = i; j < n; ++j)
		{
			if (abs(A[pivot][i]) < abs(A[j][i]))
				pivot = j;
		}

		/// if current max element in column = 0, there are no solution
		if(abs(A[pivot][i]) < EPS)
			return {Matrix (),Matrix ()};

		/// swapping matrices
		A[i].swap(A[pivot]);
		L[i].swap(L[pivot]);

		for (int j = i + 1; j < n; ++j)
		{
			double mul = A[j][i] / A[i][i];
			L[j][i] = mul;
			for (int k = 0; k < n; ++k)
			{
				A[j][k] -= (A[i][k] * mul);
			}
		}
	}

	/** fills diagonal elements with one **/
	for (int i = 0; i < n; ++i) L[i][i] = 1;

	return {L, A};
}
#endif //NUMERICALMETHODSLABS_LAB1_H











