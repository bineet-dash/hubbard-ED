#ifndef _EDLIB_H_INCLUDED_
#define _EDLIB_H_INCLUDED_

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include <lapacke.h>
#include "common_globals.h"

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue);
int bintoint(VectorXi v);
VectorXd seminvert(VectorXi  v);
long int choose(int x);
int periodic(int base, int addendum, int limit); //limit= limit starting the array from 1
bool diagonalize(MatrixXd Ac, VectorXd& lambdac, MatrixXd vc);
double find_free_energy(double temperature, vector<double> eigenvalues);

#endif
