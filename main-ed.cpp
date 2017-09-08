#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "edlib.h"

using namespace std;
using namespace Eigen;

int size;

class basis{
  int x; int spin;
public:
  basis(){x=spin=0;}
  basis(int b, int s){x=b; spin=s;}
  int get_x(){return x;}
  int get_spin(){return spin;}
  void attach_spin(int s){spin =s;}
  void output(void){cout << inttobin(x,size).transpose() << "\t \t" << spin << endl;}
};

void vector_out(std::vector<basis> v)
{
  for(auto it=v.begin(); it!=v.end(); it++)
    (*it).output();
}

void select_spin(std::vector<basis> master, std::vector<basis>& v, int spin)
{
  for(auto it=master.begin(); it!=master.end(); it++)
  {
    if((*it).get_spin()==spin)  v.push_back(*it);
  }
}

int annhilate( VectorXi v, int index, int sigma)
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
    if(inttobin(i,size).sum()==size)
    {
      int spin=seminvert(inttobin(i,size)).sum();
      half_filling.push_back(basis(i,spin));
    }
  }

  std::vector<basis> v_singlet;
  select_spin(half_filling, v_singlet, 0);
  //vector_out(v_singlet);

}
