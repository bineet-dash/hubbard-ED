#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
// #include "edlib.h"
// #include "common.h"
int size;
float t=1;
float U;
float epsilon= 1e-4;
using namespace std;
using namespace Eigen;

float find_free_energy(float temperature, vector<float> eigenvalues)
{
  float partition_func = 0;
  std::sort (eigenvalues.begin(), eigenvalues.end());
  float unruly_free_energy= 0;

  if(isinf(exp(-eigenvalues.at(0)/temperature)))
  {
    unruly_free_energy += eigenvalues.at(0);
    transform(eigenvalues.begin(), eigenvalues.end(), eigenvalues.begin(), bind1st(plus<float>(),-eigenvalues.at(0)));
  }
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++)
  {
    partition_func += exp(-(*it)/temperature);
  }
  cout << unruly_free_energy << " \t" << partition_func << endl;
  float free_energy = -unruly_free_energy/temperature + (partition_func);
  return free_energy;
}

float get_mu(float temperature, std::vector<float> v)
{
  sort (v.begin(), v.end());
  float bisection_up_lim = v.back();
  float bisection_low_lim = v.front();

  float mu;
  float no_of_electrons; int count=0;


  for(int i=0; i<=1; i++)
  {
    int it=0; no_of_electrons=0;  count++;
    cout << "Loop:" << count << "\n----------------------\n";
    mu = 0.5*(bisection_low_lim+bisection_up_lim) ;

    while(no_of_electrons<= size)
    {
      float fermi_func = 1/(exp((v[it]-mu)/temperature)+1);
      // std::cout << fermi_func << '\t';
      no_of_electrons += fermi_func;
      it++;
    }
    if(abs(v[it-1]-mu) < epsilon){cout << "exact " << v[it] << ", mu=" << mu << endl;  cout << no_of_electrons << endl; return mu; break;}
    else if(v[it-1]< mu-epsilon) {bisection_up_lim=mu; i=0; cout << no_of_electrons << endl;}
    else if(v[it-1]> mu+epsilon) {bisection_low_lim=mu; i=0; cout << no_of_electrons << endl; }
    cout << "\n-----------------------------\n";
  }

}

int main(int argc, char* argv[])
{
  if(argc!=3) exit(1);

  cout << "Enter size and U (as per input file): ";
  cin >> size >> U;

  std::vector<float> eigenvalues;
  ifstream fin(argv[1]);
    while (true)
    {
      float eival;
      fin >> eival;
      if( fin.eof() ) break;
      eigenvalues.push_back(eival);
    }

  std::cout << "Enter temperature: "; float temperature; cin >> temperature;
  cout << "The chemical potential is: " << get_mu(temperature, eigenvalues) << endl;

  // float initial_temp, final_temp, temperature_step;
  // cout << "Enter initial temperature, final temperature, temperature_step: ";
  // cin >> initial_temp >> final_temp >> temperature_step;
  //
  //  ofstream outfilefreeenergy(argv[2]);
  // for(float temperature = initial_temp; temperature < final_temp; temperature+= temperature_step)
  // {
  //   outfilefreeenergy << temperature << "\t" << get_mu(temperature, eigenvalues) << endl;
  //   //  outfilefreeenergy << temperature << " " << find_free_energy(temperature, eigenvalues) << endl;
  //   //  cout << find_free_energy(temperature,eigenvalues) << endl;
  // }
}
