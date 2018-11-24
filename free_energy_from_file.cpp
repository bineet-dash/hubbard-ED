#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

int size;
double t=1;
double U;

double find_free_energy(double temperature, vector<double> eigenvalues)
{
  std::sort (eigenvalues.begin(), eigenvalues.end());

  double partition_func = 0;
  double e_min = eigenvalues.front();
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++)
  {
    partition_func += exp(-(*it-e_min)/temperature);
  }
  double free_energy = e_min + (-temperature)*log(partition_func);
  return free_energy;
}

double find_internal_energy(double temperature, vector<double> eigenvalues)
{
  std::sort (eigenvalues.begin(), eigenvalues.end());

  double num = 0, denom=0 ; double e_min = eigenvalues.front();
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    num += exp(-(*it-e_min)/temperature)*(*it);
    denom += exp(-(*it-e_min)/temperature);
  }
  return num/denom;
}

int main(int argc, char* argv[])
{
  if(argc!=2) exit(1);

  cout << "Enter size and U (as per input file): ";
  cin >> size >> U;

  std::vector<double> eigenvalues;
  ifstream fin(argv[1]);
  while (true)
  {
    double eival;
    fin >> eival;
    if( fin.eof() ) break;
    eigenvalues.push_back(eival);
  }

  int initial_exp = -3;
  int final_exp = 0; //final_temp = 10*pow(10,final_exp)

  string filename,latticedata;
  latticedata = "_U="+to_string(int(U))+"_size_"+ to_string(size);
  filename="data/ED_results"+latticedata+".dat";
  ofstream outfile_data(filename);

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(int i=9; i>=1; i--)
    {
      double temperature = i*pow(10,j);
      outfile_data << temperature << " " << find_free_energy(temperature, eigenvalues)/double(size) << " " << find_internal_energy(temperature,eigenvalues)/double(size)  << endl;
    }
  }

  eigenvalues.clear();
  outfile_data.close();
  return 0;
}

/* double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;
    mu = 0.5*(bisection_low_lim+bisection_up_lim);

    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;
    }
    if(abs(no_of_electrons-size) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > size+epsilon)
    {
       if(abs(bisection_up_lim-v.front())<0.001){return mu; break;}
       else {bisection_up_lim=mu;}
    }
    else if(no_of_electrons < size-epsilon)
    {bisection_low_lim=mu;}
  }
}

double occupancy_check(double temperature, std::vector<double> v)
{
    sort (v.begin(), v.end());
    double mu = get_mu(temperature, v);
    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      cout << *it << " " << fermi_func << endl;
    }
} */