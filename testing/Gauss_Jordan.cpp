// TODO : rewrite it normally
#include "TestGenerator.h"

#define GN_FILE_IN STRING(generated_Gauss_Jordan.input)
#define GN_FILE_OUT STRING(generated_Gauss_Jordan.output)
#define PR_FILE_IN STRING(predefined_Gauss_Jordan.input)
#define PR_FILE_OUT STRING(predefined_Gauss_Jordan.output)

/** This file exist to test solution to SLAE using Gauss-Jordan in src/SLAE.cpp **/

using namespace std;
using namespace nm;
/** Creates stream without comments **/
void createTemp(fstream& stream)
{
	fstream toRead("temp",ios::out);
	string cur_line;
	while(getline(stream,cur_line))
	{
		cur_line = cur_line.substr(0,cur_line.find('#'));
		toRead << cur_line << '\n';
	}
	toRead.close();
}


/** prints in specified file from specified **/
void testGaussJordan(fstream& in, fstream& out)
{
	in.clear(),out.clear();
	int test = 1;
	int n,m;
	while(!in.eof())
	{
		in >> n >> m;
		Matrix curM(n,RowMatrix(m,0));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				in >> curM[i][j];
			}
		}
		Matrix toDet(n,RowMatrix(n,0));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				toDet[i][j] = curM[i][j];
			}
		}
		RowMatrix sol = nm::solveSLAEByGauss(curM);
		out << "#Test " << test++ << '\n';
		out << "Solution to this system of equation:\n";
		switch(_funcState)
		{
			case nm::ONE_SOLUTION:
				out << "There exist one solution.\n";
				break;
			case nm::INFINITY_SOLUTIONS:
				out << "There exist infinity solutions. One of them is.\n";
				break;

			case nm::NO_SOLUTION:
				out << "There don't exist any solution.\n";
				break;

			case nm::UNDEFINED_BEHAVIOUR:
				out << "Something went wrong.\n";
				break;

		}
		if(nm::_funcState == ONE_SOLUTION || nm::_funcState == INFINITY_SOLUTIONS)
		{
			for (double i : sol) out << i << ' ';
			out << '\n';
		}
		out << "Determinant: " << nm::getDeterminant(toDet) << '\n';
//		out << "Time: \n";
//		out << "Relative error: \n";
		out << "\n";
	}
}






