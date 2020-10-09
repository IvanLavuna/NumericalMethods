//
// Created by Ivan on 06.10.20.
//

#ifndef NUMERICALMETHODSLABS_LAB1_H
#define NUMERICALMETHODSLABS_LAB1_H

// TODO: rewrite it with some namespace if it becomes too large
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>

// TODO : make some defines to simplify
//  syntax of setting or unsetting flags

/// Some defines to work more easily with flags
#define SET_UNDEFINED_BEHAVIOUR state = UNDEFINED_BEHAVIOUR
#define SET_ONE_SOLUTION state 		  = ONE_SOLUTION
#define SET_INFINITY_BEHAVIOUR state  = INFINITY_SOLUTIONS
#define SET_NO_SOLUTION state 		  = NO_SOLUTION

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

// TODO : rewrite everything with this structures
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
enum  FLAG {
	ONE_SOLUTION = 0,
	INFINITY_SOLUTIONS ,
	UNDEFINED_BEHAVIOUR,
	NO_SOLUTION
} state; /// variable that defines state of function


/** helpful functions **/
void printMatrix(Matrix& m)
{
	for (auto & i : m)
	{
		for (double j : i)
		{
			cout <<std::fixed << j << " ";
		}
		cout << endl;
	}
}
void printRowMatrix(RowMatrix& m)
{
	for (double i : m)
	{
		cout << i << " ";
	}
	cout << '\n';
}

/**
 * Returns true if matrix is summetric otherwise returns false
 * **/
bool isSymmetric(Matrix & A)
{
	if(A.size() != A[0].size()) return false;
	for (int i = 0; i < A.size(); ++i)
	{
		for (int j = i+1; j < A[i].size(); ++j)
		{
			if(A[i][j] != A[j][i]) return false;
		}
	}
	return true;
}

/**
 * Transparent matrix A is matrix whose rows equals to columns of
 * same matrix A(swapping rows with columns)
 * **/
Matrix getTransparentMatrix(Matrix& A)
{
	Matrix B(A[0].size(),RowMatrix(A.size(),0));

	for (int i = 0; i < A.size(); ++i)
	{
		for (int j = 0; j < A[i].size(); ++j)
		{
			B[j][i] = A[i][j];
		}
	}
	return B;
}

/**
 *
 * @param A -> matrix that will be copied
 * @param n -> copy first n rows
 * @param m -> copy first m columns
 * @return matrix A[0 : n - 1][0 : m - 1]
 */
Matrix copyPart(Matrix& A,int n,int m)
{
	Matrix B(n,RowMatrix(m,0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			B[i][j] = A[i][j];
	}
	return B;
}
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
pair<int,RowMatrix> solveSLAEByGauss(Matrix matrix)
{
	int n = matrix.size();
	int m = static_cast<int>(matrix[0].size()) - 1;

	vector<double> solution(m);

	vector<bool> arb(m, false);

	for (int i = 0; i < std::min(n,m); ++i)
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

	int ind = m - 1;

	for (int i = m - 1; i >= 0 ; --i)
	{
		double x_i = matrix[i][m];
		for (int j = m - 1; j > i ; --j)
		{
			x_i -= (matrix[i][j] * solution[j]);
		}
		solution[ind--] = x_i;
	}

	/** if mistake of some answer is bigger than some
	 * EPS , than there are no any solution**/
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		for (int j = 0; j < m; ++j)
			sum += solution[j] * matrix[i][j];
		if (abs (sum - matrix[i][m]) > EPS)
			return {0,solution};
	}

	/** checking for infinitely many solutions **/
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
double getDeterminant(Matrix matrix)
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
int getRank(Matrix matrix)
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

/** returns P L and U matrices in PA = LU factorization **/
P_L_U get_P_L_U(Matrix A)
{
	if(A.size() != A[0].size())
	{
		throw std::invalid_argument("Matrix must be squared!");
	}
	int n = A.size();

	/** creates and fills matrix with zeros **/
	Matrix L(n,vector<double>(n,0));
	Matrix P(n,RowMatrix(n,0));
	for (int i = 0; i < n; ++i) P[i][i] = 1;


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
			return {Matrix (),Matrix (),Matrix ()};

		/// swapping matrices
		A[i].swap(A[pivot]);
		L[i].swap(L[pivot]);
		P[i].swap(P[pivot]);

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

	return {P, L, A};
}

/**
 * @brief
 * 		multiplies two matrices
 * @exception
 * 		returns exception if number of columns of
 * 		first matrix A is not equal to number of rows
 * 		in matrix B
 * 		time complexity is O(n^3)
 */
Matrix multiply(Matrix& A,Matrix& B)
{
	if(A[0].size() != B.size())
	{
		throw std::invalid_argument("You can't multiply such matrices!");
	}

	Matrix C(A.size(),vector<double>(B[0].size(),0));

	for (int i = 0; i < A.size(); ++i)
	{
		for (int j = 0; j < B[0].size(); ++j)
		{
			for (int k = 0; k < B.size(); ++k)
			{
				C[i][j] += (A[i][k] * B[k][j]);
			}
		}
	}
	return C;
}


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
RowMatrix solveSLAEByLUFactorization(Matrix& matrix)
{
	if(matrix.size() != matrix[0].size() - 1)
	{
		throw std::invalid_argument("Matrix A must be quadratic!");
	}
	int n = matrix.size();

	Matrix A(n,RowMatrix(n));
	Matrix B(n,RowMatrix(1,0)); /// ColumnMatrix

	/** copying matrix into A and B**/
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = matrix[i][j];

		B[i][0] = matrix[i][n];
	}


	P_L_U matrices = get_P_L_U(A);
	Matrix& P = matrices.P,&L = matrices.L, &U = matrices.U;

	/// if factorization not exist, single solution also not exist
	if(matrices.L.empty()) return RowMatrix(0);

	/// actual algorithm
	B = multiply(P,B);

	RowMatrix C(n,0); /// temporary matrix
	RowMatrix X(n,0); /// answer matrix

	for (int i = 0; i < n; ++i)
	{
		double c_i = B[i][0];
		for (int j = 0; j < i; ++j)
		{
			c_i -= (C[j] * L[i][j]);
		}
		C[i] = c_i;
	}

	for (int i = n - 1; i >= 0 ; --i)
	{
		double x_i = C[i];
		for (int j = n - 1; j > i ; --j)
		{
			x_i -=(U[i][j] * X[j]);
		}
		x_i /= U[i][i];
		X[i] = x_i;
	}

	return X;
}

/**
 *	finds inverse of a matrix using Gauss-Jordan
 *	algorithm.
 *	@returns zero matrix if there is not exist inverse matrix
 *	Gives @exception if matrix is not quadratic
 **/
Matrix getInverseMatrix(Matrix A)
{
	if(A.size() != A[0].size())
		throw std::invalid_argument("Matrix must be squared!");

	int n = A.size();

	/** Initially I is E matrix,but after algorithm it becomes Inverse **/
	Matrix I(n,RowMatrix(n,0));
	for(int i = 0; i < n; i++) I[i][i] = 1;

	for (int i = 0; i < n; ++i)
	{
		int pivot = i;
		for (int j = i; j < n; ++j)
		{
			if( abs(A[j][i]) > abs(A[pivot][i]))
				pivot = j;
		}

		if(abs(A[pivot][i]) < EPS) /// there do not exist I matrix
			return Matrix ();

		A[i].swap(A[pivot]);
		I[i].swap(I[pivot]);

		for (int j = 0; j < n; ++j)
		{
			if(j == i) continue;
			double mul = A[j][i] / A[i][i];
			for (int k = 0; k < n; ++k)
			{
				A[j][k] -= (A[i][k] * mul);
				I[j][k] -= (I[i][k] * mul);
			}

		}

	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			I[i][j] /= A[i][i];
	}
	return I;
}


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
L_U get_L_U_factorizationByCholesky(Matrix A)
{
	if(!isSymmetric(A))
		throw std::invalid_argument("Matrix must be symmetric!");

	int n = A.size();
	Matrix L(n,RowMatrix(n,0));

	/// using formula proposed by Wikipedia, I managed to get
	/// this factorization
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			if(i == j)
			{
				L[i][i] = A[i][i];
				for (int k = 0; k <= i - 1; ++k)
				{
					L[i][i] -= L[i][k] * L[i][k];
				}
				if(L[i][i] < 0)
				{
					SET_UNDEFINED_BEHAVIOUR;
					return L_U();
				}
				L[i][i] = sqrt(L[i][i]);
			}
			else
			{
				L[i][j] = A[i][j];
				for (int k = 0; k <= j - 1; ++k)
				{
					L[i][j] -= L[i][k] * L[j][k];
				}
				if(abs(L[j][j]) < EPS)
				{
				SET_NO_SOLUTION;
				return L_U();
				}
				L[i][j] /= L[j][j];
			}
		}
	}

	Matrix U = getTransparentMatrix(L);

	SET_ONE_SOLUTION;

	return {L, U};
}


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
RowMatrix solveSLAEByCholesky(Matrix& M)
{
	Matrix A = copyPart(M,M.size(),(int)M[0].size() - 1);
	if(!isSymmetric(A))
		throw std::invalid_argument("Matrix must be symmetric to be solved by"
							  " Cholesky!");

	int n = A.size();

	RowMatrix B(n,0);
	for (int i = 0; i < n; ++i) B[i] = M[i][n];

	L_U LU = get_L_U_factorizationByCholesky(A);
	if(state != ONE_SOLUTION)
	{
		SET_UNDEFINED_BEHAVIOUR;
		return RowMatrix();
	}

	Matrix& L = LU.L, &U = LU.U;
	RowMatrix C(n,0); /// temporary matrix
	RowMatrix X(n,0); /// answer matrix

	for (int i = 0; i < n; ++i)
	{
		double c_i = B[i];
		for (int j = 0; j < i; ++j)
		{
			c_i -= (C[j] * L[i][j]);
		}
		C[i] = c_i / L[i][i];
	}

	for (int i = n - 1; i >= 0 ; --i)
	{
		double x_i = C[i];
		for (int j = n - 1; j > i ; --j)
		{
			x_i -=(U[i][j] * X[j]);
		}
		x_i /= U[i][i];
		X[i] = x_i;
	}
	SET_ONE_SOLUTION;
	return X;
}



#endif //NUMERICALMETHODSLABS_LAB1_H











