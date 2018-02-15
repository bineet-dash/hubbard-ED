#ifndef _ED_LIBRARY_H_INCLUDED_
#define _ED_LIBRARY_H_INCLUDED_

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <clocale>
#include <cassert>
#include <lapacke.h>
#include <chrono>
#include <thread>
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

int bintoint(VectorXi v, int ph)
{
  int val=0;
  for(int i=0; i<v.size(); i++) val+= v(i)*pow(2,i);
  return ph*val;
}

VectorXd seminvert(VectorXi v_in)
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

double filter(double x) {if(abs(x)<1e-7) return 0.0; else return x;}
void filter(std::vector<double>& v) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i]); }
void filter(VectorXd& v) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i));}

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

int create(VectorXi v, int index, int sigma)
{
  int ph=0; int vec_size = v.size()/2;
  for(int i=0; i<index; i++) ph += v(i)+v(i+vec_size);

  if(sigma==-1) index+=vec_size;
  if(v(index)==0) { v(index)=1; return bintoint(v,ph);}
  else return 0;
}

int create(int x, int index, int sigma)
{
  int initial_ph = (x<0)?-1:1;
  VectorXi v= inttobin(abs(x));
  return initial_ph*create(v,index,sigma);
}

int annhilate(VectorXi v, int index, int sigma)
{
  int ph=0; int vec_size = v.size()/2;
  for(int i=0; i<index; i++) ph += v(i)+v(i+vec_size);

  if(sigma==-1) index+=v.size()/2;
  if(v(index)==1)
  {
    v(index)=0;
    return bintoint(v,ph);
  }
  else
    return 0;
}

int annhilate(int x, int index, int sigma)
{
  int initial_ph = (x<0)?-1:1;
  VectorXi v= inttobin(abs(x));
  return initial_ph*annhilate(v,index,sigma);
}

const char  *ptr = NULL;
const wchar_t uparrow[] = L"\u2191";
const wchar_t downarrow[] = L"\u2193";

void vis(int u, int d)
{
  setlocale(LC_ALL, "");
  if(u==1 && d==1)      std::wcout << uparrow << downarrow;
  else if(u==1 && d==0) std::wcout << uparrow << " ";
  else if(u==0 && d==1) std::wcout << downarrow << " ";
  else                  std::wcout << "-" << " ";
}

void vis_basis(int x, char newline)
{
  VectorXi v = inttobin(x);
  freopen(ptr, "w", stdout);
  for(int i=0; i<size; i++)
  {
    vis(v(i),v(i+size)); wcout << " ";
  }
  if(newline=='y') wcout << endl;
  else wcout << " ";
  freopen(ptr, "w", stdout);
}

class basis {
  int x; double spin; int phase;
public:
  basis(){x=spin=0;}
  basis(int b, double s){x=b; spin=s;}
  int get_x(){return x;}
  int get_phase(){return phase;}
  void get_arr(char c) {vis_basis(x,c);}
  float get_spin(){return spin;}
  void attach_spin(int s){spin = s;}
  void output(void){cout << inttobin(x).transpose() << "\t \t" << spin << endl;}
};

void vector_out(std::vector<basis> v) {for(auto it=v.begin(); it!=v.end(); it++) (*it).output();}

void select_spin(std::vector<basis> master, std::vector<basis>& v, double spin)
{ for(auto it=master.begin(); it!=master.end(); it++) if((*it).get_spin()==spin)  v.push_back(*it);}

void select_half_filling(std::vector<basis>& half_filling)
{
  long int i_min,i_max; i_min=i_max=0;
  for(int i=0; i<size; i++) i_min += pow(2,i);
  for(int i=size; i<2*size; i++) i_max += pow(2,i);
  for(int i=i_min; i<=i_max; i++)
  {
    if(inttobin(i).sum()==size) half_filling.push_back(basis(i,seminvert(inttobin(i)).sum()));
  }
}

void show_time(milliseconds begin_ms, milliseconds end_ms, string s)
{
  long int t = (end_ms.count()-begin_ms.count())/1000;
  if (t<=60)
  { cout <<  s  << t << " seconds." << endl; }
  else if (t>60 && t<3600)
  {
    int minutes=int(t/60); int seconds= t%minutes;
    cout << s  << minutes << " minute and " << seconds << " seconds." << endl;
  }
  else if(t>=3600)
  {
    int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
    cout << s  << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
  }
  else
  {cout << s  << t << "time. Wrong t received.\n"; }
}

void construct_Ht(MatrixXd& Ht, std::vector<basis> v_spin)
{
  milliseconds begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  Ht = MatrixXd::Zero(v_spin.size(),v_spin.size());
  for(int a=0; a<Ht.rows(); a++)
    {
      for(int b=0; b<Ht.rows(); b++)
      {
        Ht(a,b)=0;
        for(int sigma=-1; sigma<=1; sigma+=2)//sum over sigma
        {
          for(int i=0; i<size-1; i++)   //c\dagger_i c_i+1
          {
            int temp=annhilate(v_spin.at(b).get_x(),i+1,sigma);
            (v_spin.at(a).get_x()==create(temp,i,sigma))? Ht(a,b)+= -t: Ht(a,b)+=0;
            (v_spin.at(a).get_x()==-create(temp,i,sigma))? Ht(a,b)+= t: Ht(a,b)+=0;
          }

          for(int i=0; i<size-1; i++) //c\dagger_i+1 c_i
          {
            int temp=annhilate(v_spin.at(b).get_x(),i,sigma);
            (v_spin.at(a).get_x()==create(temp,i+1,sigma))? Ht(a,b)+= -t: Ht(a,b)+=0;
            (v_spin.at(a).get_x()==-create(temp,i+1,sigma))? Ht(a,b)+= t: Ht(a,b)+=0;
          }
        }
        cout << a << " " << b << "\r";
      }
    }
  milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  show_time(begin_ms,end_ms,"Ht construction took ");
}

void construct_HU(MatrixXd& HU, std::vector<basis> v_spin)
{
  HU = MatrixXd::Zero(v_spin.size(),v_spin.size());
  for(int a=0; a<HU.rows(); a++)
  {
    VectorXi basis = inttobin(v_spin.at(a).get_x());
    for(int i=0; i<size; i++) HU(a,a) += basis(i)*basis(i+size);
    HU(a,a) *= U;
  }
}

bool diagonalize(MatrixXd Ac, VectorXd& lambdac, MatrixXd& vc)
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

  dsyevd_(&Vchar, &Uchar, &N, Ac.data(), &LDA, lambdac.data(),  WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
  vc = Ac;
  return INFO==0;
}

bool diagonalize(MatrixXd Ac, std::vector<double>& v, MatrixXd& vc)
{
  VectorXd lambdac;
  bool result = diagonalize(Ac,lambdac,vc);
  v.resize(lambdac.size());
  VectorXd::Map(&v[0], lambdac.size()) = lambdac;
  return result;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;  count++;
    mu = 0.5*(bisection_low_lim+bisection_up_lim) ;

    auto it = v.begin();
    while(no_of_electrons< double(size) && it!=v.end())
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;  it++;
    }

    if(abs(no_of_electrons-size) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > size+epsilon)
      {
         if(bisection_up_lim == v.front()) break;
         else {bisection_up_lim=mu;}
      }
    else if(no_of_electrons < size-epsilon)
     {bisection_low_lim=mu;}
  }
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
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) partition_func += exp(-(*it)/temperature);

  double free_energy = unruly_free_energy - temperature*log(partition_func);
  return free_energy;
}

#endif
