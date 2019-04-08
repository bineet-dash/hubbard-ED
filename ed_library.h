#ifndef _ED_LIBRARY_H_INCLUDED_
#define _ED_LIBRARY_H_INCLUDED_

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <clocale>
#include <cassert>
#include <lapacke.h>
#include <chrono>
#include <thread>

extern int size;
extern double t;
extern double U;

int size;
double t=1;
double U;

using namespace std;
using namespace Eigen;
using namespace std::chrono;

typedef pair<int,double> pid;

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

double filter(double x, double x0) {if(abs(x)<x0) return 0.0; else return x;}
void filter(std::vector<double>& v, double x0) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i],x0); }
void filter(VectorXd& v, double x0) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i),x0);}

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

  if(sigma==-1)
  { for(int i=vec_size-1; i>index; i--) ph += v(i)+v(i+vec_size);}
  else if(sigma==1)
  {
    for(int i=vec_size-1; i>index; i--) ph += v(i);
    for(int i=v.size()-1; i>=index+vec_size; i--) ph += v(i);
  }

  if(sigma==-1)
  {
    index+=vec_size;
  } 
  if(v(index)==0) 
  {
    v(index)=1; 
    return bintoint(v,pow(-1,ph));
  }
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
  if(sigma==-1)
  { for(int i=vec_size-1; i>index; i--) ph += v(i)+v(i+vec_size);}
  else if(sigma==1)
  {
    for(int i=vec_size-1; i>index; i--) ph += v(i);
    for(int i=v.size()-1; i>=index+vec_size; i--) ph += v(i);
  }
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==1) { v(index)=0; return bintoint(v,pow(-1,ph));}
  else return 0;
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
const wchar_t rangle[] = L"\u3009";

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
  VectorXi v = inttobin(abs(x));
  if(x<0) cout << "-|"; else cout << " |";
  freopen(ptr, "w", stdout);
  for(int i=0; i<size; i++)
  {
    vis(v(i),v(i+size)); if(i!=size-1) wcout << " ";
  }
  wcout << rangle;
  if(newline=='y') wcout << endl;
  else wcout << " ";
  freopen(ptr, "w", stdout);
}

int n_i_left(int x)
{
 VectorXi v = inttobin(x);
 int sum = 0;
 for(int i=0; i<v.size()/4; i++) 
  {
    sum += v(i)+v(i+v.size()/2);
  }
 return sum;
}

class basis {
  int x; double spin; int phase;
public:
  basis(){x=spin=0;}
  basis(int b){x=b; spin=seminvert(inttobin(b)).sum(); phase=1;}
  basis ann(int index, int sigma) {int annhilated = annhilate(x,index,sigma); return basis(annhilated);}
  void attach_spin(int s){spin = s;}
  int get_x(){return x;}
  int get_phase(){return phase;}
  void get_arr(char c) {vis_basis(x*phase,c);}
  float get_spin(){return spin;}
  void output(void){get_arr('n'); cout << " " << spin << endl;}
};

void vector_out(std::vector<basis> v) {for(auto it=v.begin(); it!=v.end(); it++) (*it).get_arr('n');}

void select_spin( std::vector<basis> master, std::vector<basis>& v, double spin)
{ 
  for(auto it=master.begin(); it!=master.end(); it++) 
  {
    if((*it).get_spin()==spin)
      v.push_back(*it);
  }
}

void select_half_filling(std::vector<basis>& half_filling)
{
  long int i_min,i_max; i_min=i_max=0;
  for(int i=0; i<size; i++) i_min += pow(2,i);
  for(int i=size; i<2*size; i++) i_max += pow(2,i);
  for(int i=i_min; i<=i_max; i++)
  {
    if(inttobin(i).sum()==size) half_filling.push_back(basis(i));
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

void construct_Ht(MatrixXd& Ht, std::vector<basis> v_spin, char verbose_pref = 'n')
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
      if(verbose_pref=='y') cout << a << " " << b << "\r";
    }
  }
  milliseconds end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  if(verbose_pref=='y') show_time(begin_ms,end_ms,"Ht construction took ");
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

void select_filling(std::vector<basis>& selected_filling, int fill)
{
  long int i_min,i_max; i_min=i_max=0;
  // if(fill<=size)
  {
    for(int i=0; i<fill; i++) i_min += pow(2,i);
    for(int i=2*size-fill; i<2*size; i++) i_max += pow(2,i);
    for(int i=i_min; i<=i_max; i++)
    {
      if(inttobin(i).sum()==fill) selected_filling.push_back(basis(i));
    }
  }
  // else
  // {
  //   for(int i=0; i<size; i++) i_min += pow(2,i);
  //   for()
  // }
}

bool sort_pid (const pid& x1, const pid& x2) {return x1.second < x2.second;}

vector <pid> rescale(const vector <pid>& eigenvalues, double mu)
{
  vector <pid> rescaled_eivals;
  for(auto const& it : eigenvalues) rescaled_eivals.push_back(make_pair(it.first, it.second-mu*it.first));
  sort(rescaled_eivals.begin(), rescaled_eivals.end(), sort_pid);
  return rescaled_eivals;
}

double n_avg(double mu, double beta, const vector<pid> & eigenvalues)
{
  double num = 0.0, denom = 0.0;
  vector <pid> eivals_minus_mu = rescale(eigenvalues, mu);
  
  for(auto const& i: eivals_minus_mu)
  {
    num += exp(-beta*(i.second-eivals_minus_mu.front().second))*i.first;
    denom += exp(-beta*(i.second-eivals_minus_mu.front().second));
  } 
  return num/denom;
}


double get_mu(double T, vector <pid> & eigenvalues, double total_fill = size, double MU_TOLERANCE=1e-3)
{
  double beta = 1/T;
  sort(eigenvalues.begin(), eigenvalues.end(), sort_pid);
  double mu_low = eigenvalues.front().second;
  double mu_high = eigenvalues.back().second;
  double mu_mid; 
  int count = 0;

  while(true)
  {
    count++; if(count>10) exit(1);
    mu_mid = (mu_low + mu_high)/2.0;
    double suggested_fill = n_avg(mu_mid, beta, eigenvalues);

    if(abs(suggested_fill-total_fill)< MU_TOLERANCE)
    {
      // cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      break;
    }
    else if(suggested_fill < total_fill-MU_TOLERANCE)
    {
      // cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      mu_low = mu_mid;
    }
    else
    {
      // cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      mu_high = mu_mid;
    }
  }
  return mu_mid;
}


double get_free_energy(vector <pid> eigenvalues, double temperature, int fill=size)
{
  double beta = 1/temperature;

  double mu = get_mu(0.01, eigenvalues);
  vector <pid> eivals_minus_mu = rescale(eigenvalues, mu);
  double rem_F = 0.0;
  for(auto const& i : eivals_minus_mu)
  {
    rem_F += exp(-beta*(i.second-eivals_minus_mu.front().second));
  } 
  double F = eivals_minus_mu.front().second - temperature*log(rem_F) + mu*fill;
  return F;
}

#endif
