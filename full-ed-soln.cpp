#include <fstream>
#include "ed_library.h"

int main(int argc, char* argv[])
{
  cout << "Enter lattice size, U: ";
  cin >> size  >> U; assert(size%2==0);

  std::vector<pid> eigenvalues;

  for(int fill = 0; fill<=2*size; fill++)
  {
    vector<basis> selected_filling;
    select_filling(selected_filling, fill);

    double spin_limit = 0.5*fill;

    for(double sz= -spin_limit; sz<=spin_limit; sz++)
    {
      std::vector<basis> v_spin;
      select_spin(selected_filling, v_spin, sz);

      if(v_spin.size()==0) continue;

      MatrixXd Ht; construct_Ht(Ht,v_spin);
      MatrixXd HU; construct_HU(HU,v_spin);
      MatrixXd H=Ht+HU;

      std::vector<double> sz_eivals; 
      MatrixXd sz_eigenvectors;
      diagonalize(H, sz_eivals, sz_eigenvectors);
      filter(sz_eivals, 1e-5);

      for(auto const& sz_eival_i: sz_eivals) eigenvalues.push_back(make_pair(fill, sz_eival_i));
    }
  }

  // sort(eigenvalues.begin(),eigenvalues.end(), sort_pid);
  // double mu = get_mu(0.01, eigenvalues);
  // cout << "mu = " << mu << endl;

  // for(auto const& it : eigenvalues) cout << it.second << endl; cout << endl << endl;
  // for(auto const& it : eigenvalues) cout << it.second-mu*it.first << endl;

  // cout << endl << endl << "free energy = " << get_free_energy(eigenvalues, 0.001) << endl << endl;

  ofstream fout; string filename;
  filename = "data/ed_free_energy_L"+to_string(size)+"_U"+to_string(int(U))+".dat";
  fout.open(filename);

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
			double temperature = i*pow(10,j);
      fout << temperature << " " << get_free_energy(eigenvalues, temperature) << endl;
    }
  }
  fout.close();
  eigenvalues.clear();

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
