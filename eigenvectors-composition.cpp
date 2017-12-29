#include <fstream>
#include "edlib.h"
#include "common_globals.h"

typedef std::vector<pair<int,double>> eivec;

bool sort_eigenspectrum( pair<double,eivec> p1,  pair<double,eivec> p2) {return p1.first < p2.first;}
bool sort_eivec( pair <int, double> p1, pair <int, double> p2) {return abs(p1.second) > abs(p2.second);}

int main(int argc, char* argv[])
{
  cout << "Enter lattice size and U: ";
  cin >> size >> U; assert(size%2==0);

  cout << "Enter total spin (m_s): "; int spin; cin >> spin;
  if(abs(spin)>0.5*size) {cout << "Maximum m_s= " << 0.5*size << ". Exiting.\n"; exit(1);}

  vector<basis> half_filling; select_half_filling(half_filling);
  std::vector<basis> v_spin; select_spin(half_filling, v_spin, spin);
  cout << "Spin: " << spin << " sector\n----------------\nsize=" << v_spin.size() << endl;

  MatrixXd Ht; construct_Ht(Ht,v_spin);
  MatrixXd HU; construct_HU(HU,v_spin);
  MatrixXd H=Ht+HU;

  VectorXd ith_spin_eivals; MatrixXd ith_eigenvectors;
  diagonalize(H, ith_spin_eivals, ith_eigenvectors);
  filter(ith_spin_eivals);

  eivec ith_eigenvectors_listed;
  vector < pair<double,eivec> > eigenspectrum;

  for(int i=0; i<ith_spin_eivals.size(); i++)
  {
    for(int j=0; j<ith_spin_eivals.size(); j++) ith_eigenvectors_listed.push_back(make_pair(v_spin[j].get_x(),ith_eigenvectors(j,i)));
    eigenspectrum.push_back(make_pair(ith_spin_eivals(i),ith_eigenvectors_listed));
    ith_eigenvectors_listed.clear();
  }

  sort(eigenspectrum.begin(),eigenspectrum.end(),sort_eigenspectrum);

  ofstream fout; string filename;
  filename = "data/spin"+to_string(spin)+"_eivals_size"+to_string(size)+"_U"+to_string(int(U))+".txt"; fout.open(filename);
  for(auto it=eigenspectrum.begin(); it!=eigenspectrum.end(); it++) fout << (*it).first << endl;

  for(; ;)
  {
    cout << endl << "Enter an eigenvalue: ";
    int count=0; double d; cin >> d;
    for(auto it=eigenspectrum.begin(); it!=eigenspectrum.end(); it++)
    {
      if(abs((*it).first-d) < 0.01)
      {
        count++; cout << endl << "Eigenvector: " << count  << "\n=============================\n";
        sort((*it).second.begin(),(*it).second.end(),sort_eivec);
        for(int j=0; j<(*it).second.size(); j++)
         {
           // cout << ((*it).second)[j].first << " ";
           if(pow((((*it).second)[j].second),2)>1e-4)
           {
             vis_basis(((*it).second)[j].first,'n');
             cout << filter(pow((((*it).second)[j].second),2)) << endl;
           }
         }
      }
    }
    if(count==0) cout << "Not an eigenvalue!\n";
  }

  v_spin.clear(); eigenspectrum.clear();
  return 0;
}
