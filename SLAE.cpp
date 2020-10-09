//
// Created by ivan on 06.10.20.
//

#include "SLAE.h"


int main()
{
	Matrix matrix = {{2, 5, 7},
					{6, 3, 4},
					{5, -2, -3}};
	try
	{
		Matrix I = getInverseMatrix(matrix);
		Matrix B = multiply(I,matrix);
		printMatrix(I);
		cout << endl;
		printMatrix(B);

	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what() << '\n';
	}


}