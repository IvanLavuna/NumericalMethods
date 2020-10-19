//
// Created by ivan on 18.10.20.
//

#include "Matrix.h"
using std::cout;


Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) :
	Matrix(list.size(),list.begin()->size())
{
	int i = 0,j;
	for(auto& item : list)
	{
		j = 0;
		if(item.size() != _m)
			throw std::invalid_argument("Argument must have valid size of columns in each row!\n");
		for(auto & el : item)
		{
			(*this)[i][j] = el;
			++j;
		}
		++i;
	}
}

Matrix& Matrix::operator*=(const Matrix &M)
{
	if(this->_m != M._n)
	{
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - can not multiply such matrices"));
	}

	Matrix C(this->_n,M._m);
	for (int i = 0; i < _n; ++i)
	{
		for (int j = 0; j < M._m; ++j)
		{
			for (int k = 0; k < M._n; ++k)
				C[i][j] += ((*this)[i][k] * M[k][j]);
		}
	}
	(*this) = C;
	return (*this);
}

Matrix operator*(const Matrix &M1,const Matrix &M2)
{
	if(M1._m != M2._n)
	{
		throw std::invalid_argument(ERROR_MESSAGE("Invalid argument - can not multiply such matrices"));
	}

	Matrix C(M1._n, M2._m);
	for (int i = 0; i < M1._n; ++i)
	{
		for (int j = 0; j < M2._m; ++j)
		{
			for (int k = 0; k < M2._n; ++k)
				C[i][j] += M1[i][k] * M2[k][j];
		}
	}
	return C;
}

Matrix &Matrix::operator+=(const Matrix &M)
{
	if(M.rowSz() != this->rowSz() || M.columnSz() != this->columnSz())
	{
		throw std::invalid_argument(ERROR_MESSAGE("Dimensions of Matrix must be eual!"));
	}
	for (int i = 0; i < this->_n; ++i)
	{
		for (int j = 0; j < this->_m; ++j)
		{
			(*this)[i][j] += M[i][j];
		}
	}
	return *this;
}

Matrix &Matrix::operator-=(const Matrix &M)
{
	if(M.rowSz() != this->rowSz() || M.columnSz() != this->columnSz())
	{
		throw std::invalid_argument(ERROR_MESSAGE("Dimensions of Matrix must be eual!"));
	}
	for (int i = 0; i < this->_n; ++i)
	{
		for (int j = 0; j < this->_m; ++j)
		{
			(*this)[i][j] -= M[i][j];
		}
	}
	return *this;
}


Matrix &Matrix::operator*=(double n) noexcept
{
	for (int i = 0; i < this->_n; ++i)
	{
		for (int j = 0; j < this->_m; ++j)
		{
			(*this)[i][j] *= n;
		}
	}
	return *this;
}

Matrix operator+(const Matrix &M1,const Matrix &M2)
{
	if(M1.columnSz() != M2.columnSz() || M1.rowSz() != M2.rowSz())
	{
		throw std::invalid_argument(ERROR_MESSAGE("Dimensions in Matrices must be equal!"));
	}
	Matrix C = M1;
	for (int i = 0; i < M1.rowSz(); ++i)
	{
		for (int j = 0; j < M1.columnSz(); ++j)
		{
			C[i][j] += M2[i][j];
		}
	}
	return C;
}

Matrix operator-(const Matrix &M1,const Matrix &M2)
{
	if(M1.columnSz() != M2.columnSz() || M1.rowSz() != M2.rowSz())
	{
		throw std::invalid_argument(ERROR_MESSAGE("Dimensions in Matrices must be equal!"));
	}
	Matrix C = M1;
	for (int i = 0; i < M1.rowSz(); ++i)
	{
		for (int j = 0; j < M1.columnSz(); ++j)
		{
			C[i][j] -= M2[i][j];
		}
	}
	return C;
}

Matrix operator*(int n, const Matrix &M)
{
	Matrix F = M;
	for (int i = 0; i < M.rowSz(); ++i)
	{
		for (int j = 0; j < M.columnSz(); ++j)
		{
			F[i][j] *= n;
		}
	}
	return F;
}

Matrix operator*(const Matrix &M, int n)
{
	Matrix F = M;
	for (int i = 0; i < M.rowSz(); ++i)
	{
		for (int j = 0; j < M.columnSz(); ++j)
		{
			F[i][j] *= n;
		}
	}
	return F;
}


