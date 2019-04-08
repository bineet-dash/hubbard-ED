#include <fstream>
#include "ed_library.h"

void check_consistency()
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
  cin >> size >> U; assert(size%2==0);

  vector<basis> half_filling; select_half_filling(half_filling);
  std::vector<basis> v_spin;
  std::vector<double> eigenvalues;

  int spin_limit = int(0.5*size);

  for(int i= -spin_limit; i<=spin_limit; i++)
  {
    select_spin(half_filling, v_spin, i);
    cout << "Spin: " << i << " sector\n----------------\nsize=" << v_spin.size() << endl;

    MatrixXd Ht; construct_Ht(Ht,v_spin);
    MatrixXd HU; construct_HU(HU,v_spin);
    MatrixXd H=Ht+HU;

    VectorXd ith_spin_eivals_vectorxf; MatrixXd ith_eigenvectors;
    diagonalize(H, ith_spin_eivals_vectorxf, ith_eigenvectors);

    vector<double> ith_spin_eivals(ith_spin_eivals_vectorxf.data(), ith_spin_eivals_vectorxf.data()+ith_spin_eivals_vectorxf.size());
    eigenvalues.insert(eigenvalues.end(),ith_spin_eivals.begin(),ith_spin_eivals.end());

    v_spin.clear(); ith_spin_eivals.clear(); cout << endl;
  }

  sort(eigenvalues.begin(),eigenvalues.end());
  filter(eigenvalues);

  ofstream fout; string filename;
  filename = "data/eivals_size"+to_string(size)+"_U"+to_string(int(U))+".txt";
  fout.open(filename);
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) fout << *it << endl;

  check_consistency();
}
