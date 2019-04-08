#include <fstream>
#include "ed_library.h"

int main(int argc, char* argv[])
{
   // cout << "Enter lattice size and U: ";
   // cin >> size >> U; assert(size%2==0);
   int start_row, start_col, cluster_size;
   if(argc!=6) exit(1);
   start_row = std::strtol(argv[1],nullptr,0);
   start_col = std::strtol(argv[2],nullptr,0);
   cluster_size =  std::strtol(argv[3],nullptr,0);
   size =  std::strtol(argv[4],nullptr,0);
   U =  std::strtol(argv[5],nullptr,0);

   vector<basis> half_filling;
   std::vector<basis> v_spin;
   char verbose_pref = 'y';
   int spin = 0;

   select_half_filling(half_filling);
   select_spin(half_filling, v_spin, spin);

   string filename = "EDdata/ED_size="+to_string(size)+"_row="+to_string(start_row)+"_col="+to_string(start_col)+"_csize="+to_string(cluster_size)+".txt";
   ofstream fout(filename);

   for(int a = start_row; a < start_row+cluster_size; a++)
   {
      for(int b = start_col; b < start_col+cluster_size; b++)
      {
         double elem=0.0;     //Ht(a,b)
         for(int sigma=-1; sigma<=1; sigma+=2)//sum over sigma
         {
            for(int i=0; i<size-1; i++)   //c\dagger_i c_i+1
            {
               int temp=annhilate(v_spin.at(b).get_x(),i+1,sigma);
               (v_spin.at(a).get_x()==create(temp,i,sigma))? elem+= -t: elem+=0;
               (v_spin.at(a).get_x()==-create(temp,i,sigma))? elem+= t: elem+=0;
            }

            for(int i=0; i<size-1; i++) //c\dagger_i+1 c_i
            {
               int temp=annhilate(v_spin.at(b).get_x(),i,sigma);
               (v_spin.at(a).get_x()==create(temp,i+1,sigma))? elem+= -t: elem+=0;
               (v_spin.at(a).get_x()==-create(temp,i+1,sigma))? elem+= t: elem+=0;
            }
         }

         fout << elem << " ";   
         if(verbose_pref=='y') cout << a << " " << b << "\r";
      }
      fout << endl;
   }

   fout.close();
   v_spin.clear();
   return 0;
}