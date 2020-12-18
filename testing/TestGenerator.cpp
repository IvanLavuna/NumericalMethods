//
// Created by ivan on 18.10.20.
//

#include "TestGenerator.h"

void GenerateTests(std::fstream &stream,int size, int n,int m, double lower_bound, double upper_bound)
{
	for (int i = 0; i < size; ++i)
	{
		Matrix rM = GetRandomMatrix(n,m,lower_bound,upper_bound);
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

void GenerateIntTests(std::fstream &stream,int size,int n,int m, int lower_bound, int upper_bound)
{
	stream.seekp(0);
	for (int i = 0; i < size; ++i)
	{
		Matrix rM = GetIntRandomMatrix(n,m,lower_bound,upper_bound);
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

Matrix GetRandomMatrix(int n, int m,double lower_bound,double upper_bound)
{
	Matrix r_A(n,m);
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

Matrix GetIntRandomMatrix(int n, int m,int lower_bound,int upper_bound)
{
	Matrix r_A(n,m);
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

void TestGenerated(
		const std::function<RowMatrix(Matrix)>& solver,
		std::string fileName,
		int size,
		int n ,
		int m ,
		int lower_bound ,
		int upper_bound )
{

	/** creating new file to generate (with comments)**/
	fstream out(fileName + ".input",  std::ios::out),
			ans(fileName + ".output", std::ios::out);
	GenerateIntTests(out, size, n, m, lower_bound, upper_bound);
	out.close(); out.open(fileName + ".input", std::ios::in);
	/** **/

	fstream readyOut = RemoveComments(out);/** this file must be removed**/
	out.close();
	int I = 0;
	while(!readyOut.eof())
	{
		int N, M;
		readyOut >> N >> M;
		if(readyOut.eof()) break;

		ans << "Test #" << ++I << '\n';

		Matrix curMatrix(N, M);
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				readyOut >> curMatrix[i][j];
			}
		}
		/** my clock starting **/
		auto start = std::chrono::high_resolution_clock::now();
		RowMatrix curSol = solver(curMatrix);

		/** my clock ending **/
		auto end = std::chrono::high_resolution_clock::now();

		/** Solution: **/
		switch (nm::_funcState)
		{
			case nm::ONE_SOLUTION:
				ans << "There exist one solution.\n";
				break;
			case nm::INFINITY_SOLUTIONS://TODO : some bool shit here
				ans << "There exist one solutions.\n";
				break;
			case nm::UNDEFINED_BEHAVIOUR:
				ans << "Something went wrong. Probably your fault.\n";
				break;
			case nm::NO_SOLUTION:
				ans << "There don't exist any solution.\n";
				break;
		}
		if(nm::_funcState == nm::ONE_SOLUTION ||
		   nm::_funcState == nm::INFINITY_SOLUTIONS)
		{
			for (double i : curSol) ans << i << ' ';
			ans << '\n';
		}

		auto dur = end - start;
		auto f_secs = std::chrono::duration_cast<std::chrono::duration<float>>(dur);
		ans << "Time : " <<std::fixed <<  f_secs.count() << " seconds" << '\n';

		std::pair<double,double> error = GetError(curMatrix,curSol);
		ans << "Absolute error: " << std::setprecision(15) << error.first << '\n';
		ans << "Relative error: " << std::setprecision(15) << error.first << error.second << '\n';
		ans << '\n';
	}
	readyOut.close(),ans.close();

}

void TestPredefined(
		const std::function<RowMatrix(Matrix)>& solver,
		std::string fileName)
{
	/** creating new file to generate (with comments)**/
	fstream out(fileName + ".input",  std::ios::in),
			ans(fileName + ".output", std::ios::out);
	/** **/

	fstream readyOut = RemoveComments(out);/** this file must be removed**/
	out.close();
	int I = 0;
	while(!readyOut.eof())
	{
		int N, M;
		readyOut >> N >> M;
		if(readyOut.eof()) break;

		ans << "Test #" << ++I << '\n';

		Matrix curMatrix(N, M);
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				readyOut >> curMatrix[i][j];
			}
		}
		/** my clock starting **/
		auto start = std::chrono::high_resolution_clock::now();
		RowMatrix curSol = solver(curMatrix);

		/** my clock ending **/
		auto end = std::chrono::high_resolution_clock::now();

		/** Solution: **/
		switch (nm::_funcState)
		{
			case nm::ONE_SOLUTION:
				ans << "There exist one solution.\n";
				break;
			case nm::INFINITY_SOLUTIONS://TODO : some bool shit here
				ans << "There exist one solutions.\n";
				break;
			case nm::UNDEFINED_BEHAVIOUR:
				ans << "Something went wrong. Probably your fault.\n";
				break;
			case nm::NO_SOLUTION:
				ans << "There don't exist any solution.\n";
				break;
		}
		if(nm::_funcState == nm::ONE_SOLUTION ||
		   nm::_funcState == nm::INFINITY_SOLUTIONS)
		{
			for (double i : curSol) ans << i << ' ';
			ans << '\n';
		}

		auto dur = end - start;
		auto f_secs = std::chrono::duration_cast<std::chrono::duration<float>>(dur);
		ans << "Time : " <<std::fixed <<  f_secs.count() << " seconds" << '\n';

		std::pair<double,double> error = GetError(curMatrix,curSol);
		ans << "Absolute error: " << std::setprecision(15) << error.first << '\n';
		ans << "Relative error: " << std::setprecision(15) << error.first << error.second << '\n';
		ans << '\n';
	}
	readyOut.close(), ans.close();
}

fstream RemoveComments(fstream &file)
{
	fstream temp("temp", std::ios::out);
	std::string cur_line;
	while(getline(file,cur_line))
	{
		cur_line = cur_line.substr(0,cur_line.find('#'));
		temp << cur_line << '\n';
	}
	temp.close();
	temp.open("temp", std::ios::in);
	return temp;
}

/// 1 - absolute er, 2 - relative er
std::pair<double, double> GetError(Matrix &M, RowMatrix &X)
{
	RowMatrix nB(M.rowSz());
	double abs_er = 0, rel_er = 0;
	for (int i = 0; i < M.rowSz(); ++i)
	{
		for (int j = 0; j < M.columnSz() - 1; ++j)
			nB[i] -= (X[j] * M[i][j]);

		if(std::abs(nB[i] - M[i][M.columnSz() - 1]) < abs_er)
		{
			abs_er = std::abs(nB[i] - M[i][M.columnSz() - 1]);
			rel_er = std::abs(nB[i] - M[i][M.columnSz() - 1]) * 100./nB[i];
		}
	}
	return {abs_er, rel_er};

}










