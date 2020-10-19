//
// Created by ivan on 18.10.20.
//

#include "TestGenerator.h"

void generate_tests(std::fstream &stream,int size, int n,int m, double lower_bound, double upper_bound)
{
	for (int i = 0; i < size; ++i)
	{
		Matrix rM = getRandomMatrix(n,m,lower_bound,upper_bound);
		stream << "#Test " << i + 1 << '\n';
		stream << n << " " << m << '\n';
		for (int j = 0; j < n; ++j)
		{
			for (int k = 0; k < m; ++k)
			{
				stream << std::setw(10) << rM[j][k] << " ";
			}
			stream << '\n';
		}
	}
	stream << '\n';
}

void generate_int_tests(std::fstream &stream,int size,int n,int m, int lower_bound, int upper_bound)
{
	for (int i = 0; i < size; ++i)
	{
		Matrix rM = getIntRandomMatrix(n,m,lower_bound,upper_bound);
		stream << "#Test " << i + 1 << '\n';
		stream << n << " " << m << '\n';
		for (int j = 0; j < n; ++j)
		{
			for (int k = 0; k < m; ++k)
			{
				stream << std::setw(6) << rM[j][k] << " ";
			}
			stream << '\n';
		}
	}
	stream << '\n';
}

Matrix getRandomMatrix(int n, int m,double lower_bound,double upper_bound)
{
	Matrix r_A(n,RowMatrix(m,0));
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution<> uni(lower_bound,upper_bound); // guaranteed unbiased
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			r_A[i][j] = uni(rng);
		}
	}
	return r_A;
}

Matrix getIntRandomMatrix(int n, int m,int lower_bound,int upper_bound)
{
	Matrix r_A(n,RowMatrix(m,0));
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(lower_bound, upper_bound); // define the range
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			r_A[i][j] = int(distr(gen));
		}
	}
	return r_A;
}








