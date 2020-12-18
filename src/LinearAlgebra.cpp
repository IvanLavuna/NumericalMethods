//
// Created by ivan on 06.10.20.
//
#include "pch.h"
#include "LinearAlgebra.h"

namespace nm
{
	RowMatrix solveSLAEByGauss(Matrix matrix)
	{
		int n = matrix.size();
		int m = static_cast<int>(matrix[0].size()) - 1;

		RowMatrix solution(m);

		vector<bool> arb(m, false);

		for (int i = 0; i < std::min(n, m); ++i)
		{

			int pivot = i;
			/** searching for pivot element **/
			for (int j = i; j < n; ++j)
			{
				if (abs(matrix[pivot][i]) < abs(matrix[j][i]))
					pivot = j;
			}

			if (abs(matrix[pivot][i]) < EPS)
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

		for (int i = std::min(n - 1, m - 1); i >= 0; --i)
		{
			double x = matrix[i][m];
			for (int j = std::min(n - 1, m - 1); j > i; --j)
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
			if (abs(sum - matrix[i][m]) > EPS)
			{
				SET_NO_SOLUTION;
				return RowMatrix();
			}
		}

		/** checking for infinitely many solutions **/
		for (int i = 0; i < m; ++i)
		{
			if (!arb[i])
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
		if (matrix.size() != matrix[0].size())
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
			if (std::abs(matrix[pivot][i]) < EPS)
				return 0;

			if (pivot != i)
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

			if (std::abs(matrix[pivot][i]) < EPS)
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

			for (int j = i + 1; j < n; ++j)
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
		if (matrix.size() != matrix[0].size())
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
		}

		int n = matrix.size();

		/// creating matrix already filled with zeros
		Matrix L(n, n);

		/// filling diagonals with ones
		for (int i = 0; i < n; ++i) L[i][i] = 1;

		/// starting particularly like in Gauss elimination algorithm
		for (int i = 0; i < n; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				if (abs(matrix[i][i]) < EPS)
				{
					SET_NO_SOLUTION;
					return {Matrix(0,0), Matrix(0,0)};
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
		if (A.size() != A[0].size())
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
		}
		int n = A.size();

		/** creates and fills matrix with zeros**/
		Matrix L(n, n);


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
			if (abs(A[pivot][i]) < EPS)
			{
				SET_NO_SOLUTION;
				return {Matrix(0,0), Matrix(0,0)};
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
		if (A.size() != A[0].size())
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
		}
		int n = A.size();

		/** creates and fills matrix with zeros **/
		Matrix L(n,n);
		Matrix P(n, n);
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
			if (abs(A[pivot][i]) < EPS)
			{
				SET_NO_SOLUTION;
				return {Matrix(0,0), Matrix(0,0), Matrix(0,0)};
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

	RowMatrix solveSLAEByLUFactorization(Matrix &matrix)
	{
		if (matrix.size() != matrix[0].size() - 1)
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
		}
		int n = matrix.size();

		Matrix A(n, n);
		Matrix B(n, 1); /// ColumnMatrix

		/** copying matrix into A and B**/
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				A[i][j] = matrix[i][j];

			B[i][0] = matrix[i][n];
		}


		P_L_U matrices = nm::get_P_L_U(A);
		Matrix &P = matrices.P, &L = matrices.L, &U = matrices.U;

		/// if factorization not exist, single solution also not exist
		if (matrices.L.empty())
		{
			SET_NO_SOLUTION;
			return RowMatrix(0);
		}

		/// actual algorithm
		B = P * B;

		RowMatrix C(n, 0); /// temporary matrix
		RowMatrix X(n, 0); /// answer matrix

		for (int i = 0; i < n; ++i)
		{
			double c_i = B[i][0];
			for (int j = 0; j < i; ++j)
			{
				c_i -= (C[j] * L[i][j]);
			}
			C[i] = c_i;
		}

		for (int i = n - 1; i >= 0; --i)
		{
			double x_i = C[i];
			for (int j = n - 1; j > i; --j)
			{
				x_i -= (U[i][j] * X[j]);
			}
			x_i /= U[i][i];
			X[i] = x_i;
		}

		SET_ONE_SOLUTION;
		return X;
	}

	Matrix getInverseMatrix(Matrix A)
	{
		if (A.size() != A[0].size())
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be quadratic!"));
		}
		int n = A.size();

		/** Initially I is E matrix,but after algorithm it becomes Inverse **/
		Matrix I(n, n);
		for (int i = 0; i < n; i++) I[i][i] = 1;

		for (int i = 0; i < n; ++i)
		{
			int pivot = i;
			for (int j = i; j < n; ++j)
			{
				if (abs(A[j][i]) > abs(A[pivot][i]))
					pivot = j;
			}

			if (abs(A[pivot][i]) < EPS) /// there do not exist I matrix
			{
				SET_NO_SOLUTION;
				return Matrix(0,0);
			}
			A[i].swap(A[pivot]);
			I[i].swap(I[pivot]);

			for (int j = 0; j < n; ++j)
			{
				if (j == i) continue;
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
		if (!nm::isSymmetric(A))
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be symmetric!"));
		}
		int n = A.size();
		Matrix L(n, n);

		/// using formula proposed by Wikipedia, I managed to get
		/// this factorization
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j <= i; ++j)
			{
				if (i == j)
				{
					L[i][i] = A[i][i];
					for (int k = 0; k <= i - 1; ++k)
					{
						L[i][i] -= L[i][k] * L[i][k];
					}
					if (L[i][i] < 0)
					{
						SET_UNDEFINED_BEHAVIOUR;
						return {Matrix(0,0),Matrix(0,0)};
					}
					L[i][i] = sqrt(L[i][i]);
				} else
				{
					L[i][j] = A[i][j];
					for (int k = 0; k <= j - 1; ++k)
					{
						L[i][j] -= L[i][k] * L[j][k];
					}
					if (abs(L[j][j]) < EPS)
					{
						SET_NO_SOLUTION;
						return {Matrix(0,0),Matrix(0,0)};
					}
					L[i][j] /= L[j][j];
				}
			}
		}

		Matrix U = nm::TransparentMatrix(L);

		SET_ONE_SOLUTION;
		return {L, U};
	}

	RowMatrix solveSLAEByCholesky(Matrix &M)
	{
		Matrix A = nm::copyPart(M, M.size(), (int) M[0].size() - 1);
		if (!nm::isSymmetric(A))
		{
			SET_UNDEFINED_BEHAVIOUR;
			throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - Matrix must be symmetric "
													  "to be solved by Cholesky!"));
		}
		int n = A.size();

		RowMatrix B(n, 0);
		for (int i = 0; i < n; ++i) B[i] = M[i][n];

		L_U LU = nm::get_L_U_factorizationByCholesky(A);
		if (_funcState != ONE_SOLUTION)
		{
			SET_UNDEFINED_BEHAVIOUR;
			return RowMatrix();
		}

		Matrix &L = LU.L, &U = LU.U;
		RowMatrix C(n, 0); /// temporary matrix
		RowMatrix X(n, 0); /// answer matrix

		for (int i = 0; i < n; ++i)
		{
			double c_i = B[i];
			for (int j = 0; j < i; ++j)
			{
				c_i -= (C[j] * L[i][j]);
			}
			C[i] = c_i / L[i][i];
		}

		for (int i = n - 1; i >= 0; --i)
		{
			double x_i = C[i];
			for (int j = n - 1; j > i; --j)
			{
				x_i -= (U[i][j] * X[j]);
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
		int m = (int) A[0].size() - 1;

		RowMatrix X(m, 0);
		vector<bool> arb(m, false);

		for (int i = 0; i < n && i < m; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				double root = (sqrt(A[i][i] * A[i][i] + A[j][i] * A[j][i]));

				/** if elem is already = 0 we don't need to do anything **/
				if (abs(A[j][i]) < EPS) continue;
				arb[i] = true;

				double c = A[i][i] / root;
				double s = A[j][i] / root;

				for (int k = i; k <= m; ++k)
				{
					double tmp = c * A[i][k] + s * A[j][k];
					A[j][k] = -s * A[i][k] + c * A[j][k];
					A[i][k] = tmp;
				}
			}
		}

		for (int i = std::min(n - 1, m - 1); i >= 0; --i)
		{
			double x = A[i][m];
			for (int j = std::min(n - 1, m - 1); j > i; --j)
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
			if (abs(sum - A[i][m]) > EPS)
			{
				SET_NO_SOLUTION;
				return RowMatrix();
			}
		}

		// TODO: this function works wrong test

		for (int i = 0; i < m; ++i)
		{
			if (!arb[i])
			{
				SET_INFINITY_SOLUTIONS;
				return X;
			}
		}

		SET_ONE_SOLUTION;
		return X;

	}

	double conditionNumber(const Matrix &A)
	{
		Matrix inv_A = nm::getInverseMatrix(A);
		double cond = nm::MatrixNorm(A) * nm::MatrixNorm(inv_A);
		return cond;
	}

	double MatrixNorm(const Matrix &A)
	{
		double maxSum = 0;
		int n = A.size(), m = A[0].size();
		for (int i = 0; i < n; ++i)
		{
			double cur_sum = 0;
			for (int j = 0; j < m; ++j)
			{
				cur_sum += abs(A[i][j]);
			}
			maxSum = std::max(maxSum, cur_sum);
		}
		return maxSum;
	}

	double VectorNorm(const RowMatrix &A, int p)
	{
		if (p == INF)
		{
			double norm = 0;
			for (auto &item: A)
				norm = std::max(norm, abs(item));

			return norm;
		}

		double norm = 0;
		for (auto &item: A)
			norm += pow(abs(item), p);

		return pow(norm, 1. / p);
	}

	ComplexRow solveSLAEWithComplexElements(ComplexMatrix matrix)
	{
		int n = matrix.size();
		int m = static_cast<int>(matrix[0].size()) - 1;

		ComplexRow solution(m);

		vector<bool> arb(m, false);

		for (int i = 0; i < std::min(n, m); ++i)
		{

			int pivot = i;
			/** searching for pivot element **/
			for (int j = i; j < n; ++j)
			{
				if (abs(matrix[pivot][i]) < abs(matrix[j][i]))
					pivot = j;
			}

			if (abs(matrix[pivot][i]) < EPS)
				continue;

			arb[i] = true;

			matrix[i].swap(matrix[pivot]);
			std::complex<double> div = matrix[i][i];

			for (int j = i; j <= m; ++j)
			{
				matrix[i][j] /= div;
			}

			for (int j = i + 1; j < n; ++j)
			{
				std::complex<double> mul = matrix[j][i];
				for (int k = 0; k <= m; ++k)
				{
					matrix[j][k] -= (mul * matrix[i][k]);
				}
			}

		}

		for (int i = std::min(n - 1, m - 1); i >= 0; --i)
		{
			std::complex<double> x = matrix[i][m];
			for (int j = std::min(n - 1, m - 1); j > i; --j)
			{
				x -= solution[j] * matrix[i][j];
			}
			solution[i] = x;
		}

		/** if mistake of some answer is bigger than some
		 * EPS , than there are no any solution**/
		for (int i = 0; i < n; ++i)
		{
			std::complex<double> sum = 0;
			for (int j = 0; j < m; ++j)
				sum += solution[j] * matrix[i][j];
			if (abs(sum - matrix[i][m]) > EPS)
			{
				SET_NO_SOLUTION;
				return ComplexRow();
			}
		}

		/** checking for infinitely many solutions **/
		for (int i = 0; i < m; ++i)
		{
			if (!arb[i])
			{
				SET_INFINITY_SOLUTIONS;
				return solution;
			}
		}

		SET_ONE_SOLUTION;
		return solution;
	}

	RowMatrix SolveSLAEByJacobi(Matrix M)
	{
		int n = M.rowSz(), m = M.columnSz() - 1;
		RowMatrix X(m,0);
		Matrix C(n,m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				C[i][j] = -(M[i][j]/M[i][i]);
			}
		}

		double q = MatrixNorm(C);
		double p1 = (1 - q) * EPS / q;
		double p2;
		do
		{
			RowMatrix curX(m,0);

			for (int j = 0; j < m; ++j) curX[j] = M[j][m];

			for (int j = 0; j < n; ++j)
			{
				for (int k = 0; k < j; ++k)
					curX[j] -= (M[j][k] * X[k]);

				for (int k = j + 1; k < m; ++k)
					curX[j] -= (M[j][k] * X[k]);

				curX[j] /= M[j][j];
			}

			RowMatrix DX(m,0); /// delta X

			for (int i = 0; i < m; ++i)
				DX[i] = abs(curX[i] - X[i]);

			p2 = VectorNorm(DX);

			X = curX;

		} while(abs(p2) >= abs(p1));

		SET_ONE_SOLUTION;
		return X;

	}

	RowMatrix SolveSLAEByGaussSeidel(Matrix M)
	{
		int n = M.rowSz(), m = M.columnSz() - 1;
		RowMatrix X(m,0);
		Matrix C(n,m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				C[i][j] = -(M[i][j] / M[i][i]);
			}
		}

		double q = MatrixNorm(C);
		double p1 = abs((1 - q) * EPS / q);
		double p2;
		do
		{
			RowMatrix curX(m,0);

			for (int j = 0; j < m; ++j) curX[j] = M[j][m];

			for (int j = 0; j < n; ++j)
			{
				for (int k = 0; k < j; ++k)
					curX[j] -= (M[j][k] * curX[k]);

				for (int k = j + 1; k < m; ++k)
					curX[j] -= (M[j][k] * X[k]);

				curX[j] /= M[j][j];
			}

			RowMatrix DX(m,0); /// delta X

			for (int i = 0; i < m; ++i)
				DX[i] = abs(curX[i] - X[i]);

			p2 = VectorNorm(DX);

			X = curX;
		} while(abs(p2) > abs(p1));

	SET_ONE_SOLUTION;

	return X;

	}

	double GetMaxModEigenValue(const Matrix& M)
	{
		if(M.columnSz() != M.rowSz())
			throw std::invalid_argument(ERROR_MESSAGE("Matrix must be quadratic(n x n)!"));

		int n = M.size();
		double pr_e_val, cur_e_val = 0, pr_y_val = 1;
		Matrix Y(n,1); /** starting column Matrix **/
		for (int i = 0; i < n; ++i) Y[i][0] = 1;

		do
		{
			pr_e_val = cur_e_val;
			Y = M * Y;
			cur_e_val = Y[0][0] / pr_y_val;
			pr_y_val = Y[0][0];
		} while(abs(pr_e_val - cur_e_val) > EPS);

		return cur_e_val;
	}

	double GetMinEigenValue(const Matrix& M)
	{
		if(M.columnSz() != M.rowSz())
			throw std::invalid_argument(ERROR_MESSAGE("Matrix must be quadratic(n x n)!"));

		Matrix In = getInverseMatrix(M);
		double l2 = GetMaxModEigenValue(In);
		return 1. / l2;
	}

	double GetSecondMaxModEigenValue(const Matrix &M)
	{
		if(M.columnSz() != M.rowSz())
			throw std::invalid_argument(ERROR_MESSAGE("Matrix must be quadratic(n x n)!"));

		int n = M.size();
		double e_max = GetMaxModEigenValue(M);

		Matrix A = M;

		for (int j = 0; j < n; ++j)  A[j][j] -= e_max;

		return GetMaxModEigenValue(A) + e_max;
	}

	RowMatrix GetAllEigenValuesByJacobi(Matrix M)
	{
		if(!isSymmetric(M))
			throw std::invalid_argument(ERROR_MESSAGE("Matrix must be symmetric!"));

		int n = M.size();
		RowMatrix EigenValues; EigenValues.reserve(n);

		while(true)
		{
			double max_val = M[0][1];
			int I = 0, J = 1;
			double O = M_PI / 4; /// teta angle
			Matrix _M = M;

			/** finding max by mod non-diagonal element **/
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if(i == j) continue;
					if(abs(max_val) < abs(M[i][j]))
					{
						max_val = M[i][j];
						I = i, J = j;
					}
				}
			}

			if(abs(max_val) < EPS) break;

			/** calculating angle teta **/
			if(M[J][J] - M[I][I] != 0)
				O = 0.5 * atan(max_val / (M[J][J] - M[I][I]));

			double s = sin(O), c = cos(O);

			/// some complex computation just NOT to calculate Transparent(G) * M * G
			_M[I][I] = c * c * M[I][I] - 2 * s * c * max_val + s * s * M[J][J];
			_M[J][J] = s * s * M[I][I] + 2 * s * c * max_val + c * c * M[J][J];
			_M[I][J] = (c * c - s * s) * max_val + s * c * (M[I][I] - M[J][J]);
			_M[J][I] = _M[I][J];
			for (int k = 0; k < n; ++k)
			{
				if(k == I || k == J) continue;
				_M[I][k] = c * M[I][k] - s * M[J][k];
				_M[k][I] = _M[I][k];
				_M[J][k] = s * M[I][k] + c * M[J][k];
				_M[k][J] = _M[J][k];
			}

			M = _M;
		}

		for (int i = 0; i < n; ++i) EigenValues.push_back(M[i][i]);
		return EigenValues;
	}

	pair<RowMatrix, Matrix> GetAllEigenValuesAndEigenVectorsByJacobi(Matrix M)
	{
		if(!isSymmetric(M))
			throw std::invalid_argument(ERROR_MESSAGE("Matrix must be symmetric!"));

		int n = M.size();

		RowMatrix EigenValues; EigenValues.reserve(n);
		Matrix U = IdentityMatrix(n); /// Columns of this matrix will be EigenVectors

		while(true)
		{
			double max_val = M[0][1];
			int I = 0, J = 1;
			double O = M_PI / 4; /// teta angle
			Matrix G = IdentityMatrix(n); /// matrix of turns of Givens
			Matrix _M = M;

			/** finding max by mod non-diagonal element **/
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if(i == j) continue;
					if(abs(max_val) < abs(M[i][j]))
					{
						max_val = M[i][j];
						I = i, J = j;
					}
				}
			}

			if(abs(max_val) < EPS) break;

			/** calculating angle teta **/
			if(M[J][J] - M[I][I] != 0)
				O = 0.5 * atan(max_val / (M[J][J] - M[I][I]));

			double s = sin(O), c = cos(O);

			G[I][J] = s, G[J][I] = -s, G[I][I] = c, G[J][J] = c;

			/// some complex computation just NOT to calculate Transparent(G) * M * G
			_M[I][I] = c * c * M[I][I] - 2 * s * c * max_val + s * s * M[J][J];
			_M[J][J] = s * s * M[I][I] + 2 * s * c * max_val + c * c * M[J][J];
			_M[I][J] = (c * c - s * s) * max_val + s * c * (M[I][I] - M[J][J]);
			_M[J][I] = _M[I][J];
			for (int k = 0; k < n; ++k)
			{
				if(k == I || k == J) continue;
				_M[I][k] = c * M[I][k] - s * M[J][k];
				_M[k][I] = _M[I][k];
				_M[J][k] = s * M[I][k] + c * M[J][k];
				_M[k][J] = _M[J][k];
			}

			M = _M;
			U *= G;
		}

		for (int i = 0; i < n; ++i) EigenValues.push_back(M[i][i]);

		return {EigenValues, U};
	}

	double GetLagrangianValueFunction(vector<pair<double, double>>& points, double X)
	{
		SET_UNDEFINED_BEHAVIOUR;

		int n = points.size();

		double res = 0, ra, rb;

		for (int i = 0; i < n; i++)
		{
			ra = rb = 1;
			for (int j = 0; j < n; j++)
				if (i != j)
				{
					ra *= X - points[j].first;
					rb *= points[i].first - points[j].first;
				}
			res += ra * points[i].second / rb;
		}
		SET_ONE_SOLUTION;
		return res;
	}
}



