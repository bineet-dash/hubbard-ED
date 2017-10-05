#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include "common.h"

using namespace std;
using namespace Eigen;


VectorXi inttobin(int theValue)
{
  VectorXi v(2*size);
  for (int i = 0; i < 2*size; ++i)  v(i) = theValue & (1 << i) ? 1 : 0;
  return v;
}

int bintoint(VectorXi v)
{
    int val=0;
    for(int i=0; i<v.size(); i++) val+= v(i)*pow(2,i);
    return val;
}

VectorXf seminvert(VectorXi  v_in)
{
  assert(v_in.size()%2==0);
  VectorXf v = VectorXf::Zero(v_in.size());
  for(int i=0; i < v.size()/2; i++)        v_in(i)==1? v(i)= 0.5:v(i)= 0;
  for(int i=v.size()/2; i < v.size(); i++) v_in(i)==1? v(i)=-0.5:v(i)=0;
  return v;
}

long int choose(int x)
{
  long int prod =1;
  for(int i=1; i<=x; i++) prod = prod*(x+i)/i;
  return prod;
}

int periodic(int base, int addendum, int limit) //limit= limit starting the array from 1
{
  int result=base+addendum;
  if(result>=limit){result= result%(limit);}
  else if(result < 0){
    while(result < 0)
    {
      result+=limit;
    }
    result= result%(limit);
  }
  return result;
}


float find_free_energy(float temperature, vector<float> eigenvalues)
{
  float partition_func = 0;
  std::sort (eigenvalues.begin(), eigenvalues.end());
  float unruly_free_energy= 0;
  if(isinf(exp(-eigenvalues.at(0)/temperature)))
  {
    unruly_free_energy += eigenvalues.at(0);
    transform(eigenvalues.begin(), eigenvalues.end(), eigenvalues.begin(), bind1st(plus<float>(),-eigenvalues.at(0)));
  }
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++)
  {
    partition_func += exp(-(*it)/temperature);
  }
  float free_energy = unruly_free_energy - temperature*log(partition_func);
  return free_energy;
}
