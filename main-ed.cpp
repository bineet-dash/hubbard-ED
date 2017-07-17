#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "edlib.h"
using namespace std;
using namespace Eigen;

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
  cout << "Half-filling states: " << count << endl;

  int singlet_count=0;
  for(int i=0; i<no_of_combinations; i++)
     if (seminvert(inttobin(arr[i],size)).sum()==0) singlet_count++;

  cout << "Singlets: " << singlet_count << endl;

  VectorXi singlets(singlet_count); int tmp=0;
  for(int i=0; i<no_of_combinations; i++)
     {
       if (seminvert(inttobin(arr[i],size)).sum()==0)
     {
       singlets(tmp)=arr[i]; tmp++;
     }
    }


}
