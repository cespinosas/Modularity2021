#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using namespace std;


int main()
{
  int nu_exp = 600;
  
  system("mkdir -p 25_opt/redes");
  system("mkdir -p primer_opt/redes");
  system("mkdir -p pars");
  system("mkdir -p wtrack");
  system("mkdir -p dead");
  system("mkdir -p trunc");
  system("mkdir -p gen_dist");
  
  Alea jacta(3);
  Basics bas(jacta);

  ofstream fs;
  bas.open_ofstream(fs, "comandos.txt");
  
  for(int i = 0; i < nu_exp; i++)
      fs << "./chari " << i << endl;
  
  system("chmod +x comandos.txt");
  return 0;
}
