#ifndef _SPALIB_H_INCLUDED_
#define _SPALIB_H_INCLUDED_


#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue, int size);
VectorXi seminvert(VectorXi  v);
long int choose(int x);

#endif
