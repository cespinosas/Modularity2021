#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphi.h"
#include "evolvei.h"
#include "muestreoi.h"
#include "fitnessi.h"


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using namespace std;

int main()
{
    
  int seed_ances = 24125;
  int nu_exp = 600;
  
  Alea jacta(seed_ances);
  Basics bas(jacta);
  
  int nu_sfiles = 4;
  string directorio[nu_sfiles];
  directorio[1] = "primer_opt";
  directorio[2] = "dead";
  directorio[3] = "trunc";
  
  int count;
  double time1, inftrunc;
  string padead, headers;
  ifstream opf1t;
  ifstream fsd, fst;
  
  ofstream f_opf1t;
  ofstream f_fsd, f_fst;
  
  bas.open_ofstream(f_opf1t, directorio[1] + "/" + bas.inttostring(seed_ances) + "_time.txt");
    
  bas.open_ofstream(f_fsd, directorio[2] + "/" + bas.inttostring(seed_ances)  + "_dead.txt");
  
  bas.open_ofstream(f_fst, directorio[3] + "/" + bas.inttostring(seed_ances) + "_trunc.txt");
  f_fst << "red" << "\t" << "iteracion" << "\t" << "individuo" << "\t" << "generacion";
  
  for(int cl_ances = 0; cl_ances < nu_exp; cl_ances++){
    time1 = 0;
    
    bas.open_ifstream(opf1t, directorio[1] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances) + "_time.txt");
    opf1t >> time1;
    f_opf1t << time1 << endl;;
    opf1t.close();
    
    bas.open_ifstream(fsd, directorio[2] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_dead.txt");
    fsd >> padead;
    if(padead != "NA")
        f_fsd << cl_ances;
    fsd.close();
        
    bas.open_ifstream(fst, directorio[3] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_trunc.txt");
    fst >> headers >> headers >> headers;
    count = 0;
    while (fst >> inftrunc){
        if(count%3 == 0)
            f_fst << endl << cl_ances << "\t";
        if(count%3 != 2)
            f_fst << inftrunc << "\t";
        else
            f_fst << inftrunc;
        count++;
    }
    fst.close();  
  }
  
  f_opf1t.close();
  f_fsd.close();
  f_fst.close();
  
  cout << 1 << endl;
  return 0;

}
