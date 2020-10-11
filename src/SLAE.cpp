//
// Created by ivan on 06.10.20.
//

#include "SLAE.h"

using namespace nm;

void printMatrix(Matrix &m)
{
	for (auto & i : m)
	{
		for (double j : i)
			cout << std::fixed << j << " ";
		cout << '\n';
	}
}

void printRowMatrix(RowMatrix &m)
{
	for (double i : m)
		cout << i << " ";
	cout << '\n';
}

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

Matrix getTransparentMatrix(RowMatrix &A)
{
	Matrix B(A.size(),RowMatrix(1,0));

	for (int i = 0; i < A.size(); ++i)
	{
		B[i][0] = A[i];
	}
	return B;
}

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

RowMatrix solveSLAEByGauss(Matrix matrix)
{
	int n = matrix.size();
	int m = static_cast<int>(matrix[0].size()) - 1;

	RowMatrix solution(m);

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

	for (int i = std::min(n - 1,m - 1); i >= 0 ; --i)
	{
		double x = matrix[i][m];
		for (int j = std::min(n - 1,m - 1); j > i; --j)
		{
			x -= solution[j] * matrix[i][j];
		}
		solution[i] = x;
	}

	/** if mistake of some answer is bigger than some
	 * EPS , than there are no any solution**/
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		for (int j = 0; j < m; ++j)
			sum += solution[j] * matrix[i][j];
		if (abs (sum - matrix[i][m]) > EPS)
		{
			SET_NO_SOLUTION;
			return RowMatrix ();
		}
	}

	/** checking for infinitely many solutions **/
	for (int i = 0; i < m; ++i)
	{
		if(!arb[i])
		{
			SET_INFINITY_SOLUTIONS;
			return solution;
		}
	}

	SET_ONE_SOLUTION;
	return solution;

}

double getDeterminant(Matrix matrix)
{
	SET_ONE_SOLUTION;
	if(matrix.size() != matrix[0].size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Matrix must be quadratic!"));
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

int getRank(Matrix matrix)
{
	SET_ONE_SOLUTION;

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

L_U LU_factorization(Matrix matrix)
{
	if(matrix.size() != matrix[0].size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
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
				SET_NO_SOLUTION;
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
	SET_ONE_SOLUTION;
	return {L, matrix};
}

L_U PA_LU_factorization(Matrix A)
{
	if(A.size() != A[0].size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
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
		{
			SET_NO_SOLUTION;
			return {Matrix(), Matrix()};
		}
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

	SET_ONE_SOLUTION;
	return {L, A};
}

P_L_U get_P_L_U(Matrix A)
{
	if(A.size() != A[0].size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
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
		{
			SET_NO_SOLUTION;
			return {Matrix(), Matrix(), Matrix()};
		}
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

	SET_ONE_SOLUTION;
	return {P, L, A};
}

Matrix multiply(Matrix& A,Matrix& B)
{
	if(A[0].size() != B.size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - You can't multiply such matrices!"));
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
	SET_ONE_SOLUTION;
	return C;
}

RowMatrix solveSLAEByLUFactorization(Matrix &matrix)
{
	if(matrix.size() != matrix[0].size() - 1)
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
	}
	int n = matrix.size();

	Matrix A(n, RowMatrix(n));
	Matrix B(n,RowMatrix(1,0)); /// ColumnMatrix

	/** copying matrix into A and B**/
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = matrix[i][j];

		B[i][0] = matrix[i][n];
	}


	P_L_U matrices = nm::get_P_L_U(A);
	Matrix& P = matrices.P,&L = matrices.L, &U = matrices.U;

	/// if factorization not exist, single solution also not exist
	if(matrices.L.empty())
	{
		SET_NO_SOLUTION;
		return RowMatrix(0);
	}

	/// actual algorithm
	B = nm::multiply(P,B);

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

	SET_ONE_SOLUTION;
	return X;
}

Matrix getInverseMatrix(Matrix A)
{
	if(A.size() != A[0].size())
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
	}
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
		{
			SET_NO_SOLUTION;
			return Matrix();
		}
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

	SET_ONE_SOLUTION;
	return I;
}

L_U get_L_U_factorizationByCholesky(Matrix A)
{
	if(!nm::isSymmetric(A))
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be symmetric!"));
	}
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

	Matrix U = nm::getTransparentMatrix(L);

	SET_ONE_SOLUTION;
	return {L, U};
}

RowMatrix solveSLAEByCholesky(Matrix &M)
{
	Matrix A = nm::copyPart(M,M.size(),(int)M[0].size() - 1);
	if(!nm::isSymmetric(A))
	{
		SET_UNDEFINED_BEHAVIOUR;
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be symmetric "
												  "to be solved by Cholesky!"));
	}
	int n = A.size();

	RowMatrix B(n,0);
	for (int i = 0; i < n; ++i) B[i] = M[i][n];

	L_U LU = nm::get_L_U_factorizationByCholesky(A);
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

RowMatrix solveSLAEByMethodOfTurns(Matrix A)
{
	int n = A.size();
	int m = (int)A[0].size() - 1;

	RowMatrix X(m,0);
	vector<bool> arb(m, false);

	for (int i = 0; i < n && i < m; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double root = (sqrt(A[i][i] * A[i][i] + A[j][i] * A[j][i]));

			/** if elem is already = 0 we don't need to do anything **/
			if(abs(A[j][i]) < EPS) continue;
			arb[i] = true;
			double c = A[i][i] / root;
			double s = A[j][i] / root;

			for (int k = i; k <= m; ++k)
			{
				double tmp = c * A[i][k] + s * A[j][k];
				A[j][k] =  - s * A[i][k] + c * A[j][k];
				A[i][k] = tmp;
			}
		}
	}

	for (int i = std::min(n - 1,m - 1); i >= 0 ; --i)
	{
		double x = A[i][m];
		for (int j = std::min(n - 1,m - 1); j > i; --j)
		{
			x -= X[j] * A[i][j];
		}
		X[i] = x / A[i][i];
	}

	/** if mistake of some answer is bigger than some
 	*  EPS , than there are no any solution **/
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		for (int j = 0; j < m; ++j)
			sum += X[j] * A[i][j];
		if (abs (sum - A[i][m]) > EPS)
		{
			SET_NO_SOLUTION;
			return RowMatrix();
		}
	}


	for (int i = 0; i < m; ++i)
	{
		if(!arb[i])
		{
			SET_INFINITY_SOLUTIONS;
			return X;
		}
	}

	SET_ONE_SOLUTION;
	return X;

}
