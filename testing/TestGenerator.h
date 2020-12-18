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
#include "Matrix.h"
#include <fstream>
#include <random>
#include <iomanip>
#include <functional>
#include <chrono>
using std::fstream;
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
void GenerateTests(std::fstream &stream,int size, int n,int m, double lower_bound, double upper_bound);

/** generates size number of test which are integers with range [from:to] **/
void GenerateIntTests(std::fstream &stream,int size,int n,int m, int lower_bound, int upper_bound);

/** generates and returns Matrix n x m with random values **/
Matrix GetRandomMatrix(int n, int m,double lower_bound,double upper_bound);

/** generates and returns Matrix n x m with random int values **/
Matrix GetIntRandomMatrix(int n, int m,int lower_bound,int upper_bound);

/** 1) Creates file *.input filled with random data that are coefficients of Matrix(Generates tests),
 *  2) then it inputs all generated tests and uses function solver to solve given equation,
 *  3) and outputs it to *.output file in following format
 *  `
 * There exist one/zero/infinity solutions to this problem.
 * Solution : 		`rowMatrix`
 * Time : 			`time`
 * Absolute error : 'absolute error'
 * Relative error : 'relative error'
 * `
 * Uses integer coefficients by default
 * **/
void TestGenerated(const std::function<RowMatrix(Matrix)>& solver,
				   std::string fileName,
				   int size 		 = 10,
				   int n 			 = 10,
				   int m 			 = 10, /** one more than n because we include column b**/
				   int lower_bound   = -30,
				   int upper_bound   = 30);

void TestPredefined(const std::function<RowMatrix(Matrix)>& solver,
					std::string fileName);

/** returns file called temp without comments **/
fstream RemoveComments(fstream& out);

/** return absolute ans relative errors in pair**/
std::pair<double,double> GetError(Matrix& M, RowMatrix& X);
#endif //NUMERICALMETHODSLABS_TESTGENERATOR_H










