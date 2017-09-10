  #include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "edlib.h"
#include "common.h"

using namespace std;
using namespace Eigen;

int size;
float t=1;

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


int annhilate(VectorXi v, int index, int sigma)
{
  if(sigma==-1) index+=v.size()/2;
  if(v(index)==1)
  {
    v(index)=0;
    return bintoint(v);
  }
  else
  {
    return 0;
  }
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
  {
    return 0;
  }
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
  {
    return 0;
  }
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
  {
    return 0;
  }
}

void vector_out(std::vector<basis> v)
{
  for(auto it=v.begin(); it!=v.end(); it++)
    (*it).output();
}

void select_spin(std::vector<basis> master, std::vector<basis>& v, float spin)
{
  for(auto it=master.begin(); it!=master.end(); it++)
  {
    if((*it).get_spin()==spin)  v.push_back(*it);
  }
}

int main()
{

  cout << "Enter lattice size: ";
  cin >> size;

  long int no_of_combinations,i_min,i_max;
  i_min=i_max=0;
  //no_of_combinations = choose(size);

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

  std::vector<basis> v_singlet;
  select_spin(half_filling, v_singlet, 0.);
  cout << "Singlet basis are: \n";
  vector_out(v_singlet);

  MatrixXf Ht(v_singlet.size(),v_singlet.size());

  for(int a=0; a<Ht.rows(); a++)
  {
    for(int b=0; b<Ht.rows(); b++)
    {
      Ht(a,b)=0;
      for(int sigma=-1; sigma<=1; sigma+=2)
      {
        for(int i=0; i<size; i++)
        {
           int temp=annhilate(v_singlet.at(b).get_x(),periodic(i,1,size),sigma);
           (v_singlet.at(a).get_x()==create(temp,i,sigma))? Ht(a,b)+= -t: Ht(a,b)+=0;
        }
      }
    }
  }

  cout << "Ht matrix is: \n\n" << Ht << endl;

}
