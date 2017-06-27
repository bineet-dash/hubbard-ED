#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue, int size)
{
  VectorXi v1(2*size);
  for (int i = 0; i < 2*size; ++i)  v1(i) = theValue & (1 << i) ? 1 : 0;
  return v1;
}

VectorXi seminvert(VectorXi  v)
{
  if(v.size()%2==1) exit(1);
  for(int i=v.size()/2; i < v.size(); i++) v(i)==1? v(i)=-1:v(i) =0;
  return v;
}

long int choose(int x)
{
  long int prod =1;
  for(int i=x+1; i<=2*x; i++) prod*=i;
  for(int i=1; i<=x; i++) prod /= i;
  return prod;
}


int main()
{
  cout << "Enter lattice size: ";
  int size;
  cin >> size;

  long int no_of_combinations,i_min,i_max;
  i_min=i_max=0; no_of_combinations = choose(size);


  for(int i=0; i<size; i++) i_min += pow(2,i);
  for(int i=size; i<2*size; i++) i_max += pow(2,i);

  int count = 0;
  int* arr= new (nothrow) int[no_of_combinations];

  for(int i=i_min; i<=i_max; i++)
    if(inttobin(i,size).sum()==size) {arr[count] = i; count++;}

  assert(count==no_of_combinations);
  cout << "Half-filling states: " << count << endl << endl;

  int singlet_count=0;
  for(int i=0; i<no_of_combinations; i++)
     if (seminvert(inttobin(arr[i],size)).sum()==0) singlet_count++;

  cout << "Singlets: " << singlet_count << endl;
  delete []arr;

}
