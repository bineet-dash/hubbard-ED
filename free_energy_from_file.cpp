#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
// #include "edlib.h"
#include "common.h"

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
  float free_energy = unruly_free_energy - temperature*log(partition_func);
  return free_energy;
}


int main(int argc, char* argv[])
{
  if(argc!=3) exit(1);

  std::vector<float> eigenvalues;

  ifstream fin(argv[1]);
    while (true)
    {
      float eival;
      fin >> eival;
      if( fin.eof() ) break;
      eigenvalues.push_back(eival);
    }

  float initial_temp, final_temp, temperature_step;
  cout << "Enter initial temperature, final temperature, temperature_step: ";
  cin >> initial_temp >> final_temp >> temperature_step;

  // cout << find_free_energy(initial_temp, eigenvalues) << endl;

  // ofstream outfilefreeenergy(argv[2]);
  for(float temperature = initial_temp; temperature < final_temp; temperature+= temperature_step)
  {
    // outfilefreeenergy << temperature << " " << find_free_energy(temperature, eigenvalues) << endl;
    cout << find_free_energy(temperature,eigenvalues) << endl;
  }
}
