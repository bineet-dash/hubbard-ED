#ifndef _EDLIB_H_INCLUDED_
#define _EDLIB_H_INCLUDED_


#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue, int size);
int bintoint(VectorXi v);


VectorXf seminvert(VectorXi  v);
long int choose(int x);

#endif
