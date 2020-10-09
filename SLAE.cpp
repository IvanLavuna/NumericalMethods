//
// Created by ivan on 06.10.20.
//

#include "SLAE.h"


int main()
{
	Matrix matrix = {{1,  1, 3,  3,  12},
					 {1,  5, 3,  2,  32},
					 {3,  3, 19, 11,  23},
					 {3,  2, 11,  14,  4}};
	try
	{
		RowMatrix A = solveSLAEByCholesky(matrix);
		RowMatrix B = solveSLAEByGauss(matrix).second;
		printRowMatrix(A);
		cout << endl;
		printRowMatrix(B);

	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what() << '\n';
	}


}