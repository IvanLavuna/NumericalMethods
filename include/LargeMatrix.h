//
// Created by ivan on 20.10.20.
//

#ifndef NUMERICALMETHODSLABS_LARGEMATRIX_H
#define NUMERICALMETHODSLABS_LARGEMATRIX_H

#include "Matrix.h"

/** Matrix that can represent sparse matrices
 * This data structure will only store non-zero elements
 * **/

class LargeMatrix
{
private:
	vector<double> _vArr;
	vector<int>    _lj;
	vector<int>    _li;

public:
	LargeMatrix(vector<vector<double>> & M);
	LargeMatrix(Matrix & M);
	
};
















#endif //NUMERICALMETHODSLABS_LARGEMATRIX_H
