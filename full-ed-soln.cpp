#include <fstream>
#include "ed_library.h"
#include "common_globals.h"

typedef pair<int,double> pid;

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

bool sort_pid(const pid& x1, const pid& x2){return x1.second < x2.second;}

double n_avg(double mu, double beta, vector<pid> eigenvalues)
{
  double num = 0.0, denom = 0.0;
  for(auto const& i: eigenvalues)
  {
    num += exp(-beta*(i.second-mu*i.first))*i.first;
    denom += exp(-beta*(i.second-mu*i.first));
  } 
  // cout << num << " " << denom << endl;
  return num/denom;
}

double get_mu(double T, vector <pid> eigenvalues, double total_fill = size, double MU_TOLERANCE=1e-3)
{
  double beta = 1/T;
  sort(eigenvalues.begin(), eigenvalues.end(), sort_pid);
  double mu_low = eigenvalues.front().second;
  double mu_high = eigenvalues.back().second;
  double mu_mid; 

  while(true)
  {
    mu_mid = (mu_low + mu_high)/2.0;
    double suggested_fill = n_avg(mu_mid, beta, eigenvalues);

    if(abs(suggested_fill-total_fill)< MU_TOLERANCE)
    {
      cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      break;
    }
    else if(suggested_fill < total_fill-MU_TOLERANCE)
    {
      cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      mu_low = mu_mid;
    }
    else
    {
      cout << suggested_fill << " " << mu_high << " " << mu_mid << " " << mu_low << " " << endl; 
      mu_high = mu_mid;
    }
  }
  return mu_mid;
}

vector <pid> rescale(const vector <pid>& eigenvalues, double mu)
{
  vector <pid> rescaled_eivals;
  for(auto const& it : eigenvalues) rescaled_eivals.push_back(make_pair(it.first, it.second-mu*it.first));
  sort(rescaled_eivals.begin(), rescaled_eivals.end(), sort_pid);
  return rescaled_eivals;
}

double get_free_energy( const vector <pid> & eigenvalues, double temperature, int fill=size)
{
  double beta = 1/temperature;
  sort(eigenvalues.begin(),eigenvalues.end(), sort_pid);
  double mu = get_mu(0.01, eigenvalues);
  vector <pid> eivals_minus_mu = rescale(eigenvalues, mu);
  double rem_F = 0.0;
  for(auto const& i : eivals_minus_mu) rem_F += exp(beta*(i.second-eivals_minus_mu.front().second));
  double F = eivals_minus_mu.front().second - temperature*log(rem_F);
  return F + mu*fill;
}

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
      cout << "Spin: " << sz << " sector\n----------------\nsize=" << v_spin.size() << endl;

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

  sort(eigenvalues.begin(),eigenvalues.end(), sort_pid);
  double mu = get_mu(0.01, eigenvalues);
  cout << "mu = " << mu << endl;

  for(auto const& it : eigenvalues) cout << it.second << endl; cout << endl << endl;
  for(auto const& it : eigenvalues) cout << it.second-mu*it.first << endl;

  cout << get_free_energy(eigenvalues, 0.01) << endl << endl;

  ofstream fout; string filename;
  filename = "data/full_ed_eivals_size"+to_string(size)+"_U"+to_string(int(U))+".txt";
  fout.open(filename);
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++) fout << (*it).first << " " << (*it).second << endl;
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
