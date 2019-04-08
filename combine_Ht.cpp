#include <iostream>
#include <locale>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <experimental/filesystem>
#include <Eigen/Dense>
#include <fstream>
#include "ed_library.h"

using namespace std;
using namespace Eigen;

std::string extract_ints(std::ctype_base::mask category, std::string str, std::ctype<char> const& facet)
{
  using std::strlen;

  char const *begin = &str.front(),
              *end   = &str.back();

  auto res = facet.scan_is(category, begin, end);

  begin = &res[0];
  end   = &res[strlen(res)];

  return std::string(begin, end);
}

std::string extract_ints(std::string str)
{
    return extract_ints(std::ctype_base::digit, str,
         std::use_facet<std::ctype<char>>(std::locale("")));
}

std::vector<std::string> explode(std::string const & s, char delim)
{
    std::vector<std::string> result;
    std::istringstream iss(s);

    for (std::string token; std::getline(iss, token, delim); )
    {
        result.push_back(std::move(token));
    }

    return result;
}

void get_params(string filename, int& size, int& start_row, int& start_col, int& cluster_size)
{
  int params [5]; 
  auto v = explode(filename, '_');
  for(int i=0; i<v.size(); i++)
  {
    string s = v.at(i);
    stringstream ss(extract_ints(s)); 
    ss >> params[i];
  } 
  size = params[1];
  start_row = params[2];
  start_col = params[3];
  cluster_size = params[4];
}

namespace fs = std::experimental::filesystem;

int main(int argc, char* argv[])
{
  if(argc!=2) {cout << "give matrix size\n"; exit(1);}
  int m_size = stoi(argv[1]);
  MatrixXd Ht = MatrixXd::Constant(m_size,m_size,9);

  for (auto & p : fs::directory_iterator("EDdata/"))
  {
    std::ostringstream oss;
    oss << p;
    std::string path = oss.str();
    if ( path.front() == '"' ) 
    {
      path.erase( 0, 1 ); // erase the first character
      path.erase( path.size() - 1 ); // erase the last character
    }
    cout << path << endl;
    ifstream fin(path);
    if(fin.fail())
    {
      cout << path << " not avaliable." << endl;
      continue; 
    }

    int size, start_row, start_col, csize;
    get_params(path, size, start_row, start_col, csize);

    for(int i=start_row; i<start_row+csize; i++)
    {
      for(int j=start_col; j<start_col+csize; j++)
       {
         fin >> Ht(i,j); 
        //  double a; fin >> a; cout << a; 
       }
    }
    fin.close()
  };

  // VectorXd lambda; MatrixXd vc;
  // diagonalize(Ht, lambda, vc);
  // ofstream fout("check_eivals.txt");
  // fout << lambda << endl;
  // fout.close();
  
  ofstream fout;
  fout.open("check_Ht.txt");
  fout << Ht << endl;
  fout.close();
  
  return 0;
}