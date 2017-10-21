#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include <fstream>
#include <lapacke.h>
#include <cmath>
#include <clocale>
#include "edlib.h"
#include "common.h"

using namespace std;
using namespace Eigen;

int size;
double t=1;
double U;

const char  *ptr = NULL;
const wchar_t up[] = L"\u2191";
const wchar_t down[] = L"\u2193";

void vis(int u, int d)
{
  setlocale(LC_ALL, "");
  if(u==1 && d==1)     std::wcout << up << down;
  else if(u==1&& d==0) std::wcout << up;
  else if(u==0&& d==1) std::wcout << down;
  else                 std::wcout <<"O";
}

class basis {
  int x; double spin;
public:
  basis(){x=spin=0;}
  basis(int b, double s){x=b; spin=s;}
  int get_x(){return x;}
  void get_arr(char);
  float  get_spin(){return spin;}
  void attach_spin(int s){spin =s;}
  void output(void){cout << inttobin(x).transpose() << "\t \t" << spin << endl;}
};

void basis::get_arr(char newline)
{
  VectorXi v = inttobin(x);
  freopen(ptr, "w", stdout);
  for(int i=0; i<size; i++)
  {
    vis(v(i),v(i+size)); wcout << " ";
  }
  if(newline=='y') wcout << endl;
  freopen(ptr, "w", stdout);
}

double filter(double x) {if(abs(x)<1e-4) return 0.0; else return x;}
void filter(std::vector<double>& v) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i]); }
VectorXd filter(VectorXd v) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i)); return v;}

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


int annhilate(VectorXi v, int index, int sigma)
{
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==1)
  {
    v(index)=0;
    return bintoint(v);
  }
  else
    return 0;
}

int annhilate(int x, int index, int sigma)
{
  VectorXi v=inttobin(x);
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==1)
  {
    v(index)=0;
    return bintoint(v);
  }
  else
    return 0;
}


int create(VectorXi v, int index, int sigma)
{
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==0) { v(index)=1; return bintoint(v);}
  else return 0;
}

int create(int x, int index, int sigma)
{
  VectorXi v= inttobin(x);
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==0){ v(index)=1; return bintoint(v);}
  else return 0;
}

void vector_out(std::vector<basis> v) {for(auto it=v.begin(); it!=v.end(); it++) (*it).output();}

void select_spin(std::vector<basis> master, std::vector<basis>& v, double spin)
{ for(auto it=master.begin(); it!=master.end(); it++) if((*it).get_spin()==spin)  v.push_back(*it);}

void check_consistency(double t, double U)
{
  if(size!=2) return;
  cout << "-------------------------------------\n";
  cout << "Theoretical result for eigenvalues: \n";
  cout << U/2-sqrt(U*U/4+4*t*t) << endl << 0 << " (Triplet) " << endl << U << endl << U/2+sqrt(U*U/4+4*t*t)<< endl;
}

void check_tb_validity(void)
{
  cout << "The eigenvalues obtained by theoretical TB model: \n";
  for(int i=0; i<size; i++) cout << -2*t*cos(2*M_PI*i/double(size)) << ", ";
  cout << endl;
}

int main(int argc, char* argv[])
{
  cout << "Enter lattice size and U: ";
  cin >> size >> U;
  assert(size%2==0);

  long int no_of_combinations,i_min,i_max;
  i_min=i_max=0;

  for(int i=0; i<size; i++) i_min += pow(2,i);
  for(int i=size; i<2*size; i++) i_max += pow(2,i);

  int count = 0;
  vector<basis> half_filling;

  for(int i=i_min; i<=i_max; i++)
  {
    if(inttobin(i).sum()==size)
    {
      double spin=seminvert(inttobin(i)).sum();
      half_filling.push_back(basis(i,spin));
    }
  }

  std::vector<basis> v_spin;
  std::vector<double> eigenvalues;

  int spin_limit = int(0.5*size);

  for(int i= -spin_limit; i<=spin_limit; i++)
  {
      select_spin(half_filling, v_spin, i);

      cout << "Spin: " << i << " sector\nsize=" << v_spin.size() << endl;

      MatrixXd Ht(v_spin.size(),v_spin.size());
      for(int a=0; a<Ht.rows(); a++)
        {
          for(int b=0; b<Ht.rows(); b++)
          {
            Ht(a,b)=0;
            for(int sigma=-1; sigma<=1; sigma+=2)
            {
              for(int i=0; i<size; i++)
              {
                 int temp=annhilate(v_spin.at(b).get_x(),periodic(i,1,size),sigma);
                 (v_spin.at(a).get_x()==create(temp,i,sigma))? Ht(a,b)+= -t: Ht(a,b)+=0;
              }
            }
          }
        }

      MatrixXd HU= MatrixXd::Zero(v_spin.size(),v_spin.size());
      for(int a=0; a<HU.rows(); a++)
        {
          VectorXi basis = inttobin(v_spin.at(a).get_x());
          for(int i=0; i<size; i++) HU(a,a) += basis(i)*basis(i+size);
          HU(a,a) *= U;
        }

      MatrixXd H=Ht+HU; VectorXd ith_spin_eivals_vectorxf; VectorXd ith_eigenvectors;
      diagonalize(H, ith_spin_eivals_vectorxf, ith_eigenvectors);

      vector<double> ith_spin_eivals(ith_spin_eivals_vectorxf.data(), ith_spin_eivals_vectorxf.data()+ith_spin_eivals_vectorxf.size());
      eigenvalues.insert(eigenvalues.end(),ith_spin_eivals.begin(),ith_spin_eivals.end());

      v_spin.clear();
      ith_spin_eivals.clear();
  }

    sort(eigenvalues.begin(),eigenvalues.end());
    filter(eigenvalues);

    ofstream fout; string filename;
    filename = "data/eivals_size"+to_string(size)+"_U"+to_string(U)+".txt";
    fout.open(filename);
    for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) fout << *it << endl;

    // check_consistency(1,U);
    // check_tb_validity();

}
