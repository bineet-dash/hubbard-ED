#ifndef _EDLIB_H_INCLUDED_
#define _EDLIB_H_INCLUDED_


#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include "common.h"

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue);
int bintoint(VectorXi v);
VectorXf seminvert(VectorXi  v);
long int choose(int x);
int periodic(int base, int addendum, int limit); //limit= limit starting the array from 1


#endif
