#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue)
{
  VectorXi v1(12);
  for (int i = 0; i < 12; ++i)
  {
      v1(i) = theValue & (1 << i) ? 1 : 0;
  }
  return v1;
}

int factorial(int x)
{
  int fac;
  x>0 ? fac = x*factorial(x-1) : fac = 1;
  return fac;
}

void progress_percent(float initial_i, float final_i, float i)
{
   int progress_percent;
   progress_percent=int((i)*100/(final_i-initial_i));
   cout.flush();
   cout << "\r [ "<< progress_percent+1 << "% ] ";
   for(int i=0; i<progress_percent; i++)
   {
     cout << "#";
   }
}


int main()
{
  int size;
  cin >> size;

  int no_of_combinations,i1,i2;
  no_of_combinations = i1=i2=0;

  for(int i=0; i<size; i++)
  {
    i1+= pow(2,i);
  }

  for(int i=size; i<2*size; i++)
  {
    i2+= pow(2,i);
  }

  int one_count; int count = 0;
  no_of_combinations = factorial(2*size)/pow(factorial(size),2);
  int* arr= new int[no_of_combinations];

  for(int i=i1; i<=i2; i++)
  {
    one_count = 0;
    for(int j=0; j<2*size; j++)
      if( i & (1 << j) ) one_count++;
    if(one_count==size) {arr[count] = i; count++; }
  }
  cout << count << endl << endl;

  int t1, t2;
  int singlet_count=0;
  for(int i=0; i<no_of_combinations; i++)
    {
      t1 = t2=0;
      for(int j=0; j<size; j++)
        if(arr[i] & (1 << j)) t1++;
      for(int j=size; j<2*size; j++)
        if(arr[i] & (1 << j)) t2++;
      if(t1==t2) singlet_count++ ;

    }

  cout << "Singlets: " << singlet_count << endl;

}
