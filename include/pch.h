//
// Created by ivan on 19.10.20.
//

#ifndef NUMERICALMETHODSLABS_PCH_H
#define NUMERICALMETHODSLABS_PCH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <initializer_list>


/// Some defines to work more easily with flags
#define SET_UNDEFINED_BEHAVIOUR _funcState = UNDEFINED_BEHAVIOUR
#define SET_ONE_SOLUTION _funcState		  = ONE_SOLUTION
#define SET_INFINITY_SOLUTIONS _funcState  = INFINITY_SOLUTIONS
#define SET_NO_SOLUTION _funcState 		  = NO_SOLUTION

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define ERROR_MESSAGE(DESC) ("error: " AT "\n" "Description: " DESC "\n")





#endif //NUMERICALMETHODSLABS_PCH_H
