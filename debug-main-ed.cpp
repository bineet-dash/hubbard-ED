#define EIGEN_USE_LAPACKE

#include <iostream>
#include <Eigen/Dense>
#include <cassert>
#include <fstream>
#include <lapacke.h>
#include <cmath>
#include "edlib.h"
#include "common.h"
#include <clocale>

const wchar_t up[] = L"\u2191";
const wchar_t down[] = L"\u2193";


using namespace std;
using namespace Eigen;

int size;
float t=1;
float U;

class basis {
  int x; float spin;
public:
  basis(){x=spin=0;}
  basis(int b, float s){x=b; spin=s;}
  int get_x(){return x;}
  float  get_spin(){return spin;}
  void attach_spin(int s){spin =s;}
  void output(void){cout << inttobin(x).transpose() << "\t \t" << spin << endl;}
};

void vis(int u, int d)
{
  setlocale(LC_ALL, "");
  if(u==1 && d==1)     std::wcout << up << down;
  else if(u==1&& d==0) std::wcout << up;
  else if(u==0&& d==1) std::wcout << down;
  else                 std::wcout <<"O";
}

void show_arr(int theValue)
{
  VectorXi v = inttobin(theValue);
  for(int i=0; i<size; i++)
  {
    vis(v(i),v(i+size)); wcout << " ";
  }
}

float filter(float x) {if(abs(x)<1e-4) return 0.0; else return x;}
void filter(std::vector<float>& v) {for(int i=0; i<v.size(); i++)  v[i]=filter(v[i]); }
VectorXf filter(VectorXf v) {for(int i=0; i<v.size(); i++)  v(i)=filter(v(i)); return v;}

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
  if(v(index)==0)
  {
    v(index)=1;
    return bintoint(v);
  }
  else
    return 0;
}

int create(int x, int index, int sigma)
{
  VectorXi v= inttobin(x);
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==0)
  {
    v(index)=1;
    return bintoint(v);
  }
  else
    return 0;
}

void vector_out(std::vector<basis> v)
{
  for(auto it=v.begin(); it!=v.end(); it++)
    (*it).output();
}

void select_spin(std::vector<basis> master, std::vector<basis>& v, float spin)
{
  for(auto it=master.begin(); it!=master.end(); it++)
    if((*it).get_spin()==spin)  v.push_back(*it);
}

void check_consistency(float t, float U)
{
  if(size!=2) return;
  cout << "-------------------------------------\n";
  cout << "Theoretical result for eigenvalues: \n";
  cout << U/2-sqrt(U*U/4+4*t*t) << endl << 0 << " (Triplet) " << endl << U << endl << U/2+sqrt(U*U/4+4*t*t)<< endl;
}

void check_tb_validity(void)
{
  cout << "The eigenvalues obtained by theoretical TB model: \n";
  for(int i=0; i<size; i++) cout << -2*t*cos(2*M_PI*i/float(size)) << ", ";
  cout << endl;
}

int main(int argc, char* argv[])
{
  //assert(argc>1);
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
      float spin=seminvert(inttobin(i)).sum();
      half_filling.push_back(basis(i,spin));
    }
  }

  std::vector<basis> v_spin;
  std::vector<float> eigenvalues;

  int spin_limit = int(0.5*size);
  int i =0;

      select_spin(half_filling, v_spin, i);

      MatrixXf Ht(v_spin.size(),v_spin.size());
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

      // EigenSolver <MatrixXf> es;
      // es.compute(H);

      cout << "Spin:" << i << " sector Hamiltonian: \n";
      for(int i=0; i<v_spin.size(); i++) cout << v_spin[i].get_x() << " ";
      cout << endl;
      for(int i=0; i<v_spin.size(); i++) cout << "--";
      cout << endl;
      cout << Ht << endl << endl;

      while(1==1)
      {
        const char  *ptr = NULL;

        int row, col;
        cout << "Enter row and column: ";
        cin >> row >> col;

        cout << "row array for row=" << row << ": ";
        freopen(ptr, "w", stdout);
        show_arr(v_spin.at(row).get_x());

        freopen(ptr, "w", stdout);
        wcout << "\ncol array for col=" << col << ": ";

        freopen(ptr, "w", stdout);
        show_arr(v_spin.at(col).get_x());
        wcout << endl;
        freopen(ptr, "w", stdout);
      }

      exit(1);

    //   VectorXf ith_spin_eivals_vectorxf = es.eigenvalues().real();
    //   std::vector<float> ith_spin_eivals(ith_spin_eivals_vectorxf.data(), ith_spin_eivals_vectorxf.data()+ith_spin_eivals_vectorxf.size());
    //   eigenvalues.insert(eigenvalues.end(),ith_spin_eivals.begin(),ith_spin_eivals.end());
    //
    //   {cout << "The singlet eivals are: \n" << filter(ith_spin_eivals_vectorxf).transpose() << endl << endl;}
    //
    //   v_spin.clear();
    //   ith_spin_eivals.clear();
    //
    //
    // sort(eigenvalues.begin(),eigenvalues.end());
    // filter(eigenvalues);

    // ofstream fout;
    // fout.open(argv[1]);
    // for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) cout << *it << endl;
    // check_consistency(1,U);
    //check_tb_validity();

}
