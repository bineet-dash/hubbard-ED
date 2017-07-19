#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "edlib.h"

using namespace std;
using namespace Eigen;

int main()
{
  int t; cin >> t;
  cout << inttobin(t,10).transpose() << endl;

  VectorXi v=inttobin(t,10);
  cout << bintoint(v) << endl;
}
