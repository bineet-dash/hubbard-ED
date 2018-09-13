#include <fstream>
#include <cstdlib>
#include "ed_library.h"
#include "common_globals.h"

typedef std::vector<pair<int,double>> eivec;

bool sort_eigenspectrum( pair<double,eivec> p1,  pair<double,eivec> p2) {return p1.first < p2.first;}
bool sort_eivec( pair <int, double> p1, pair <int, double> p2) {return abs(p1.second) > abs(p2.second);}

int main(int argc, char* argv[])
{
  int Nr; //Size of Ht
  cout << "Enter lattice size and Ht size: ";
  cin >> size >> Nr; assert(size%2==0);

  int spin = 0;  
  vector<basis> half_filling; select_half_filling(half_filling);
  std::vector<basis> v_spin; select_spin(half_filling, v_spin, spin);

  MatrixXd Ht = MatrixXd::Constant(Nr,Nr,15000);
  ifstream Htin("check_Ht.txt");
  for(int i=0; i<Nr; i++)
  {
    for(int j=0; j<Nr; j++)
    {
      Htin >> Ht(i,j);
    }
  }

  ofstream dataout("deln_vs_U_size_"+to_string(size)+".txt");

  for(double i=0; i<20; i+=1)
  {
    ::U=i;
    MatrixXd HU; construct_HU(HU,v_spin);
    MatrixXd H=Ht+HU;
    // MatrixXd H = Ht;

    VectorXd ith_spin_eivals; MatrixXd ith_eigenvectors;

    cout << "diagonalization started" << endl;
    diagonalize(H, ith_spin_eivals, ith_eigenvectors);
    filter(ith_spin_eivals,1e-5);
    cout << "diagonalization finished" << endl;

    eivec ith_eigenvectors_listed;
    vector < pair<double,eivec> > eigenspectrum;

    for(int i=0; i<ith_spin_eivals.size(); i++)
    {
      for(int j=0; j<ith_spin_eivals.size(); j++) ith_eigenvectors_listed.push_back(make_pair(v_spin[j].get_x(),ith_eigenvectors(j,i)));
      eigenspectrum.push_back(make_pair(ith_spin_eivals(i),ith_eigenvectors_listed));
      ith_eigenvectors_listed.clear();
    }
    sort(eigenspectrum.begin(),eigenspectrum.end(),sort_eigenspectrum);
    
    double temperature = 1e-5;
    double n_av = 0.0;
    double n2_av = 0.0;
    double e0 = eigenspectrum.at(0).first;
    double part_func = 0.0;

    for(auto it=eigenspectrum.begin(); it!=eigenspectrum.end(); it++)
    {
      double n_i_l = 0.0;
      double n_i_l2 = 0.0;
      double adjusted_eival = (*it).first-e0; 
      for(int j=0; j<(*it).second.size(); j++)
      {
          int x = ((*it).second)[j].first;
          n_i_l += n_i_left(x)*norm(((*it).second)[j].second);
          n_i_l2 += pow(n_i_left(x),2)*norm(((*it).second)[j].second);
      }
      part_func += exp(-adjusted_eival/temperature);
      n_av += exp(-adjusted_eival/temperature)*n_i_l;
      n2_av += exp(-adjusted_eival/temperature)*n_i_l2;
    }

    n_av = n_av/part_func;
    n2_av = n2_av/part_func;
    cout << "U= " << U << "\n===============\n";
    cout << n_av << " " << n2_av << endl; 
    cout << endl << n2_av - pow(n_av,2) << " " << n2_av/pow(n_av,2)-1 << endl << endl;
    dataout << U << " " << n2_av - pow(n_av,2) << " " << n2_av/pow(n_av,2)-1 << endl;
    eigenspectrum.clear();
  }

  v_spin.clear(); 
  dataout.close();
  return 0;
}

