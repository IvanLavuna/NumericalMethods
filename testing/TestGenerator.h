//
// Created by ivan on 18.10.20.
//

#ifndef NUMERICALMETHODSLABS_TESTGENERATOR_H
#define NUMERICALMETHODSLABS_TESTGENERATOR_H
//
// Created by ivan on 09.10.20.
//

/**
 * This code provides basic functions to test code in Numerical Methods Labs
 *
 * Note:
 * All values in Matrices has type DOUBLE. When i wrote function
 * generate_int_test or someting like that I actually mean that
 * values are integers with type 1.00000000;
 *
 */
#include "LinearAlgebra.h"
#include <fstream>
#include <random>
#include <iomanip>
using std::fstream;
using nm::Matrix;
using nm::RowMatrix;
#define STRING(s) #s

#define PATH STRING(/home/ivan/CLionProjects/NumericalMethodsLabs/testing/)

/**
 * generates [size] number of test with range [from:to] and
 * puts all those tests into specified file [stream]
 * in format:
 *
 * 		#Test [n]
 * 		n m
 * 		[Matrix ]
 * 		empty
 * 	n - row size of matrices
 * 	m - column size of matrices
 * 	size - number of test cases
 * **/
void generate_tests(std::fstream &stream,int size, int n,int m, double lower_bound, double upper_bound);

/** generates size number of test which are integers with range [from:to] **/
void generate_int_tests(std::fstream &stream,int size,int n,int m, int lower_bound, int upper_bound);

/** generates and returns Matrix n x m with random values **/
Matrix getRandomMatrix(int n, int m,double lower_bound,double upper_bound);

/** generates and returns Matrix n x m with random int values **/
Matrix getIntRandomMatrix(int n, int m,int lower_bound,int upper_bound);

#endif //NUMERICALMETHODSLABS_TESTGENERATOR_H










