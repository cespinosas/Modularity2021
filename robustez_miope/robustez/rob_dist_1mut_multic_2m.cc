#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphi.h"
#include "cap_evoli.h"
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
  
  int seed = 24125;
  int nodes = 24;
  int modules = 4;
  int sample_nt = 5000;
  int at_step = 500;
  int nu_exp = 5000;
  int nic = pow(2, (modules-1));

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI red(jacta);
  Cap_evolI chw(jacta);
  FitnessI lafi(jacta);

  set<set<int> > predpar;
  set<int>* i_par;
  i_par = new set<int>[modules];
  
  for(int j = 0; j < modules; j++){
    for (int i = j*(nodes/modules); i < (j+1)*(nodes/modules); i++)
        i_par[j].insert(i);
    predpar.insert(i_par[j]);
    i_par[j].clear();
  }

  
  int *s;
  s = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       if ((i%2)==0)
            s[i] = 1;
       else
            s[i] = -1;
  }
  
  int nuic;
  nuic = pow(2, modules - 1);
  int **ats;
  bas.create_array(ats, nuic, nodes);
  chw.start_all_ics(ats, nuic, modules, nodes);
  
  int *mgood, *mbad;
  int **oc, **ot;
  
  mgood = new int[nodes];
  mbad = new int[nodes];
  oc = new int*[nodes];
  ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      oc[i] = new int[modules];
      ot[i] = new int[modules];
  }

  double rob, robprom, dist;;
  
  double rob_gan, rob_per, dist_gan, dist_per, mrob_gan, mrob_per;
  int cuangan, cuanper, trunc;
  
  string directorio = "muestrario";
  
  ifstream muestra;
  ifstream inventario;
  ofstream orob;
  
  bas.open_ofstream(orob, "rob/2mic_"  + bas.inttostring(seed) + "_allmut_2m.txt");
  orob << "Rob" << "\t" << "Rob_by_dist" << "\t" << "Rob_prom" << "\t" << "Rob_gan" << "\t" << "Rob_per" << "\t" << "Rob_by_dist_gan" << "\t" << "Rob_by_dist_per" << "\t" << "Rob_prom_gan" << "\t" << "Rob_prom_per" << "\t" << "Cuan_gan" << "\t" << "Cuan_per" << "\t" << "Trunc" << endl;
  //Obtiene los genotipos extremos
  for (int i = 0; i < nu_exp; i++){
//       cout << i << endl;
      bas.open_ifstream(muestra, directorio + "/" + bas.inttostring(seed) + "_muestra_2m.txt");
      bas.open_ifstream(inventario, directorio + "/" + bas.inttostring(seed) + "_inventario_2m.txt");

      red.get_dir_nw_from_file_bignw(nodes, muestra, inventario, i, sample_nt);
      
      rob = 0;
      dist = 0;
      if(red.rob_by_dist_1mutation_cap(ats, nic, s, mgood, mbad, oc, ot, nodes, modules, predpar, at_step, rob, dist, robprom, rob_gan, rob_per, dist_gan, dist_per, mrob_gan, mrob_per, cuangan, cuanper, trunc)){
          orob << rob << "\t" << dist << "\t" << robprom << "\t" << rob_gan  << "\t" << rob_per << "\t" << dist_gan << "\t" << dist_per << "\t" << mrob_gan << "\t" << mrob_per << "\t" << cuangan << "\t" << cuanper << "\t" << trunc << endl;
      }
      else{
          cout << "[Error]: Network did not produce the modular attractors.\n";
          exit(1);
      }
      red.clear();
      
      muestra.close();
      inventario.close();
  }
  orob.close();
  
  jacta.close_rng();
  cout << 1 << endl;  
}
