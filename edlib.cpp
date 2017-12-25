#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include <lapacke.h>
#include <chrono>
#include "common_globals.h"

int size;
double t=1;
double U;

using namespace std;
using namespace Eigen;
using namespace std::chrono;

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

VectorXd seminvert(VectorXi  v_in)
{
  assert(v_in.size()%2==0);
  VectorXd v = VectorXd::Zero(v_in.size());
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

bool diagonalize(MatrixXd Ac, VectorXd& lambdac, MatrixXd vc)
{
  int N;
  if(Ac.cols()==Ac.rows())  N = Ac.cols(); else return false;

  lambdac.resize(N);
  vc.resize(N,N);

  int LDA = N;
  int INFO = 0;
  char Uchar = 'U';
  char Vchar = 'V';
  char Nchar = 'N';

  int LWORK = 5*(2*LDA*LDA+6*LDA+1);
  int LIWORK = 5*(3+5*LDA);

  VectorXd WORK(LWORK);
  VectorXi IWORK(IWORK);

  dsyevd_(&Nchar, &Uchar, &N, Ac.data(), &LDA, lambdac.data(),  WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  vc = Ac;
  return INFO==0;
}

double find_free_energy(double temperature, vector<double> eigenvalues)
{
  double partition_func = 0;
  std::sort (eigenvalues.begin(), eigenvalues.end());
  double unruly_free_energy= 0;
  if(isinf(exp(-eigenvalues.at(0)/temperature)))
  {
    unruly_free_energy += eigenvalues.at(0);
    transform(eigenvalues.begin(), eigenvalues.end(), eigenvalues.begin(), bind1st(plus<double>(),-eigenvalues.at(0)));
  }
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++)
  {
    partition_func += exp(-(*it)/temperature);
  }
  double free_energy = unruly_free_energy - temperature*log(partition_func);
  return free_energy;
}

void show_time(milliseconds begin_ms, milliseconds end_ms,string s)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout <<  s << " took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << s << " took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << s << " took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << s << " took " << t << "time. Wrong t received.\n"; }
}
