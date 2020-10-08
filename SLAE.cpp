//
// Created by ivan on 06.10.20.
//

#include "SLAE.h"


int main()
{
	Matrix matrix = { {2,1,5},
					 {4,4,-4},
					 {1,3,1}
				        };
	try
	{
		pair<Matrix ,Matrix > p = PA_LU_factorization(matrix);

		for (int i = 0; i < matrix.size(); ++i)
		{
			for (int j = 0; j < matrix.size(); ++j)
			{
				cout << p.first[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;

		for (int i = 0; i < matrix.size(); ++i)
		{
			for (int j = 0; j < matrix.size(); ++j)
			{
				cout << p.second[i][j] << " ";
			}
			cout << endl;
		}
	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what();
	}

}