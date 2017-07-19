#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <cassert>

using namespace std;
using namespace Eigen;

VectorXi inttobin(int theValue, int size)
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

// float matrix_element_ht(MatrixXf M, int ket, int bra, float t)
// {
//   int size = M.rows()/2;
//   VectorXi bin_ket= inttobin(ket,size);
//   VectorXi bin_bra= inttobin(bra,size);
//   for(int annhiliation_index=0; annhiliation_index<size; annhiliation_index++)
//   {
//     creation_index=periodic(annhiliation_index,1,size);
//
//   }
// }
//
//
// void construct_ht(MatrixXf M, float t)
// {
//   int size = M.rows()/2;
//   int creation_index, annihilation_index;
//
//   for(int sigma=-1; sigma<=1; sigma+=2)
//   {
//     for(int i=0; i<size; i++)
//      {
//        annihilation_index=i;
//        creation_index=periodic(i,1,lattice_size);
//
//        for(int ket=0; x<2*size; x++)    //matrix for (c_j\dagger c_i)
//        {
//          for(int bra=0; y<2*size; y++)
//          {
//            M(x,y)+= -t*matrixelement_ht(x,y,lattice_size,creation_index,annihilation_index,sigma);
//          }
//        }
//      }
//    }
// }

int matrixelement_h0(int x, int y, int lattice_size, int i, int j, int sigma)
{
  int* arr = new int[lattice_size]; int tmp;
  arr=int_to_arr(y, lattice_size);
  if(arr[j]==sigma && arr[i]==0)
  {
    arr[i]=sigma;
    arr[j]=0;
    tmp= arr_to_int(arr,lattice_size);
    if(x==tmp) return 1;
    else return 0;
  }
  else return 0;
  delete arr;
}



void obsolete_construct_h0(MatrixXf &M, int lattice_size, float t)
{
  int creation_index, annihilation_index;
  //cout << "here func\n";
  for(int sigma=-1; sigma<=1; sigma+=2)
  {
    //cout << "here sigma"<< sigma << " \n";
    for(int i=0; i<lattice_size; i++)
     {
       //cout << "here index" << i <<" \n";
       annihilation_index=i;
       creation_index=periodic(i,1,lattice_size);

       for(int x=0; x<2*lattice_size; x++)    //matrix for (c_j\dagger c_i)
       {
         for(int y=0; y<2*lattice_size; y++)
         {
           //cout << "matrixelement_h0 " << x << " " << y << '\n';
           M(x,y)+= -t*matrixelement_h0(x,y,lattice_size,creation_index,annihilation_index,sigma);
         }
       }
       //cout << "index " << i << " complete \n";
     }
        //cout << "sigma " << sigma << "complete\n";
   }

   //cout << "constructed!\n";

  // for(int x=0; x<2*lattice_size; x++) //matrix for (c_i\dagger c_j)
  //  {
  //    for(int y=0; y<2*lattice_size; y++)
  //    {
  //        if(M(x,y)!=0) M(y,x)=M(x,y);
  //    }
  //  }
}
