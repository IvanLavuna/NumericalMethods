#include <iostream>
#include "LinearAlgebra.h"
#include <iomanip>
using namespace nm;
using namespace std;

Matrix A {
		{4,6,0,0,0,  0},
		{3,13,2,0,   0},
		{0,3,8,5,0,  0},
		{0,0,2,7,3,  0},
		{0,0,0,5,9, -4},
		{0,0,0,0,4, 13}
};

int main()
{
	try
	{
		P_L_U b = get_P_L_U(A);

		printMatrix(b.L);
		cout << endl;
		printMatrix(b.U);
	}
	catch (std::exception& exception)
	{
		std::cerr << exception.what();
	}
}









