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

  if(isinf(exp(-eigenvalues.at(0)/temperature)))
  {
    transform(eigenvalues.begin(), eigenvalues.end(), eigenvalues.begin(), bind1st(plus<double>(),-eigenvalues.at(0)));
  }

  double partition_func = 0;
  for(auto it=eigenvalues.begin(); it!=eigenvalues.end(); it++)
  {
    partition_func += exp(-(*it)/temperature);
  }
  double free_energy = eigenvalues.at(0) + (-temperature)*log(partition_func);
  return free_energy;
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

  // for( ; ; )
  // {
  //   double temperature; cin >> temperature;
  //   double free_energy = find_free_energy(temperature, eigenvalues);
  //   cout << "free_energy = " << free_energy << endl << endl;
  // }

  int initial_exp = -3;
  int final_exp = 2;

  ofstream outfile_freeenergy; string filename,latticedata;
  latticedata = "_U"+to_string(int(U))+"_size"+ to_string(size);
  filename="data/ED_free_energy_vs_temp"+latticedata+"_T_"+to_string(int(10*pow(10,final_exp)))+".txt"; outfile_freeenergy.open(filename);

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(int i=10; i>=2; i--)
    {
      double temperature = i*pow(10,j);
      outfile_freeenergy << temperature << " " << find_free_energy(temperature, eigenvalues) << endl;
    }
  }

  eigenvalues.clear();
  outfile_freeenergy.close();
  return 0;
}

// double debug_get_mu(double temperature, std::vector<double> v)
// {
//   sort (v.begin(), v.end());
//   double bisection_up_lim = v.back();
//   double bisection_low_lim = v.front();
//
//   double mu, no_of_electrons; int count=0;
//   double epsilon = 0.000001;
//
//   for(; ;)
//   {
//     no_of_electrons=0;  count++; //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//   //  if(count>10) exit(1);
//     cout << "Loop:" << count << "\n----------------------\n";
//     mu = 0.5*(bisection_low_lim+bisection_up_lim) ;
//
//     cout << "new mu = " << mu << endl;
//
//     auto it = v.begin();
//     while(no_of_electrons< double(size) && it!=v.end())
//     {
//       double fermi_func = 1/(exp((*it-mu)/temperature)+1);
//       std::cout << fermi_func << " ";
//       no_of_electrons += fermi_func;  it++;
//     }
//
//     cout << endl;
//
//     if(abs(no_of_electrons-size) < epsilon)
//     {
//       cout << "exact " << *it << ", mu=" << mu << " " << ", no_of_electrons= " << no_of_electrons << endl;
//       return mu; break;
//     }
//     else if(no_of_electrons > size+epsilon)
//       {
//          if(bisection_up_lim == v.front()) break;
//          else {cout << "upper " << mu << " " << " " <<  no_of_electrons << endl; bisection_up_lim=mu; std::cerr << "new low =" << bisection_low_lim << " new up =" << bisection_up_lim << '\n';}
//       }
//     else if(no_of_electrons < size-epsilon)
//      { cout << "lower " << mu << " " << " " << no_of_electrons << endl; bisection_low_lim=mu;}
//     cout << "\n-----------------------------\n";
//   }
// }


// double maxval(std::vector<double> v) {sort(v.begin(),v.end()); return v.back();}
// double minval(std::vector<double> v) {sort(v.begin(),v.end()); return v.front();}
//
// double get_mu_s(double fill, std::vector<double> evl_s, double T)
// {
//   double f0, f, fl2, fr, mr, ml, rtmp, m_d;
//
//   mr = maxval(evl_s);
//   fr=0.0;
//   for(int i=0; i<evl_s.size(); i++)
//     fr += (1.0/(exp((evl_s[i]-mr)/T)+1.0));
//
//   ml= minval(evl_s);
//   fl2=0.0;
//   for(int i=0; i<evl_s.size(); i++)
//     fl2 += (1.0/(exp((evl_s[i]-ml)/T)+1.0));
//
//   m_d= 0.5*(ml+mr);
//   f=0.0;
//   for(int i=0; i<evl_s.size(); i++)
//     f += (1.0/(exp((evl_s[i]-m_d)/T)+1.0));
//
//     std::cerr << "===============" << endl;
//     cerr << " mr= " << mr  << " ml= " << ml << " m_d= " << m_d << " \n";
//     std::cerr << "================" << '\n';
//
//   int count =0;
//
//   while(abs(f-fill)>1e-2)
//   {
//     cerr << ++count << ") ";
//
//     m_d=0.5*(ml+mr);
//
//     cerr << " mr= " << mr  << " ml= " << ml << " m_d= " << m_d << " ";
//
//     f=0.0;
//     for(int i=0; i<evl_s.size(); i++)
//       f += (1.0/(exp((evl_s[i]-m_d)/T)+1.0));
//     if(f>fill)
//     {
//       mr= m_d; fr=f;
//       cerr << " /upper/ ";
//     }
//     else
//     {
//       ml=m_d; fr = f;
//       cerr << " /lower/ ";
//     }
//
//     cerr << " f= " << f << endl;
//     if(count==100) break;
//   }
//
//   return m_d;
// }

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