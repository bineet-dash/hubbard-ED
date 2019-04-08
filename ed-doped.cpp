#include <fstream>
#include "ed_library.h"
// #include <matplotlibcpp.h>

// namespace plt=matplotlibcpp;

void check_consistency();

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


int main(int argc, char* argv[])
{
  int fill;
  cout << "Enter lattice size, filling, U: ";
  cin >> size >> fill >> U; assert(size%2==0);

  vector<basis> selected_filling;
  select_filling(selected_filling, fill);
  
  std::vector<basis> v_spin;
  std::vector<double> eigenvalues;

  double spin_limit = 0.5*fill;

  for(double i= -spin_limit; i<=spin_limit; i++)
  {
    select_spin(selected_filling, v_spin, i);

    cout << "Spin: " << i << " sector\n----------------\nsize=" << v_spin.size() << endl;

    if(v_spin.size()==0) continue;

    MatrixXd Ht; construct_Ht(Ht,v_spin);
    MatrixXd HU; construct_HU(HU,v_spin);
    MatrixXd H=Ht+HU;

    // cout << Ht << endl << endl;

    std::vector<double> ith_spin_eivals; MatrixXd ith_eigenvectors;
    diagonalize(H, ith_spin_eivals, ith_eigenvectors);

    eigenvalues.insert(eigenvalues.end(),ith_spin_eivals.begin(),ith_spin_eivals.end());
    v_spin.clear(); ith_spin_eivals.clear(); cout << endl;
  }

  sort(eigenvalues.begin(),eigenvalues.end());
  filter(eigenvalues,1e-5);

  ofstream fout; string filename;
  filename = "data/eivals_size"+to_string(size)+"_filling"+to_string(fill)+"_U"+to_string(int(U))+".txt";
  fout.open(filename);
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) fout << *it << endl;
  fout.close();

  // plt::hist(eigenvalues,20);
  // plt::show();

  eigenvalues.clear();

  // check_consistency();

  return 0;
}

void check_consistency(void)
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
