//
// Created by ivan on 16.10.20.
//
#include "pch.h"
#include "Basic.h"

namespace nm
{
	FLAG _funcState = UNDEFINED_BEHAVIOUR;


	void printMatrix(Matrix &m)
	{
		for (auto &i : m)
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

	void printComplexRowMatrix(ComplexRow &m)
	{
		for (auto &i : m)
			cout << "(" << i.real() << "," << i.imag() << ")";
		cout << '\n';
	}

	bool isSymmetric(const Matrix &A)
	{
		if (A.size() != A[0].size()) return false;
		for (int i = 0; i < A.size(); ++i)
		{
			for (int j = i + 1; j < A[i].size(); ++j)
			{
				if (A[i][j] != A[j][i]) return false;
			}
		}
		return true;
	}

	Matrix TransparentMatrix(Matrix &A)
	{
		Matrix B(A.columnSz(),A.rowSz());

		for (int i = 0; i < A.size(); ++i)
		{
			for (int j = 0; j < A[i].size(); ++j)
			{
				B[j][i] = A[i][j];
			}
		}
		return B;
	}

	Matrix TransparentMatrix(RowMatrix &A)
	{
		Matrix B(A.size(), 1);

		for (int i = 0; i < A.size(); ++i)
		{
			B[i][0] = A[i];
		}
		return B;
	}

	Matrix copyPart(Matrix &A, int n, int m)
	{
		Matrix B(n, m);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				B[i][j] = A[i][j];
		}
		return B;
	}

	Matrix IdentityMatrix(int n)
	{
		Matrix A(n,n);
		for (int i = 0; i < n; ++i)  A[i][i] = 1;
		return A;
	}

}