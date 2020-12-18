//
// Created by ivan on 18.10.20.
//

#ifndef NUMERICALMETHODSLABS_MATRIX_H
#define NUMERICALMETHODSLABS_MATRIX_H

// TODO: make this class be more like "matrix" - NOT like some piece of shit
//  Also add RowMatrix and ColumnMatrix

/** Matrix class that uses int his basis vector of vectors
 *
 * **/

using std::vector;

class Matrix : public vector<vector<double>>
{
private:
	int _n; /// size of rows
	int _m; /// size of columns

public:
	Matrix(int n,int m):
		vector<vector<double>>(n,vector<double>(m,0)),
		_n(n),
		_m(m)
	{}
	explicit Matrix(std::initializer_list<std::initializer_list<double>> list);

	int rowSz()    const{ return _n; }
	int columnSz() const{ return _m; }

	/// operators
	Matrix& operator*=(const Matrix& M);
	Matrix& operator+=(const Matrix& M);
	Matrix& operator-=(const Matrix& M);
	Matrix& operator*=(double n) noexcept;

	friend Matrix operator*(const Matrix& M1, const Matrix&M2);
	friend Matrix operator+(const Matrix& M1, const  Matrix& M2);
	friend Matrix operator-(const Matrix& M1, const  Matrix& M2);
	friend Matrix operator*(const Matrix& M, int n);
	friend Matrix operator*(int n,const Matrix& M);

	friend std::ostream& operator<<(std::ostream& out,Matrix M)
	{
		for (auto& i : M)
		{
			for (auto& j : i)
			{
				out << j << ' ';
			}
			out << '\n';
		}
		return out;
	}

};


#endif //NUMERICALMETHODSLABS_MATRIX_H












