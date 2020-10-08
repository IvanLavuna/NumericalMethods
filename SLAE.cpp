//
// Created by ivan on 06.10.20.
//

#include "SLAE.h"


int main()
{
	Matrix matrix = {
					 {1,3,1,2,  6},
					 {2,1,5,12,  5},
					 {4,4,-4,13, 0},
					 {11,3,1,129,  6},
					 {4,4,-4,1203, 0}};
	try
	{
		pair<int, RowMatrix> m2 = solveSLAEByGauss(matrix);
		cout << m2.first << "\n";
		printRowMatrix(m2.second);
		cout << endl;

	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what() << '\n';
	}


}