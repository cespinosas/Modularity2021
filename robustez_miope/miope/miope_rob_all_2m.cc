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
  int samplesize = 5000;
  int num_exp = 400;
  int nodes = 24;
  int modules = 4;
  int seed = 24125;
  int conectividad = 3*nodes;
  
  int max_mut = 2*nodes;
  double wol = double(conectividad)/double(nodes*nodes);
  
  Alea jacta(seed);
  Basics bas(jacta);
  Cap_evolI cap(jacta);
  GraphI red(jacta);
  FitnessI fifu(jacta);
  
  set<set<int> > predpar;
  set<int>* i_par;
  i_par = new set<int>[modules];
  
  for(int j = 0; j < modules; j++){
    for (int i = j*(nodes/modules); i < (j+1)*(nodes/modules); i++)
        i_par[j].insert(i);
    predpar.insert(i_par[j]);
    i_par[j].clear();
  }
  delete [] i_par;

  int *s;
  s = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       if ((i%2)==0)
            s[i] = 1;
       else
            s[i] = -1;
  }
  
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
  
  double **matmod;
  bas.create_array(matmod, nodes, nodes);
  double rob;
  double last_rob, new_rob;
  
  ofstream ocon;
  bas.open_ofstream(ocon, "miope_rob/" +  bas.inttostring(seed) + "_con_2m.txt");
  ocon << "Initial" << "\t" << "Final" << endl;
    
  ofstream omod;
  bas.open_ofstream(omod, "miope_rob/" + bas.inttostring(seed) + "_mod_2m.txt");
  omod << "Initial" << "\t" << "Final" << endl;
  
  ofstream orob;
  bas.open_ofstream(orob, "miope_rob/" + bas.inttostring(seed) + "_rob_2m.txt");
  orob << "Initial" << "\t" << "Final" << endl;
  
  ofstream otime;
  bas.open_ofstream(otime, "miope_rob/" + bas.inttostring(seed) + "_time_2m.txt");
  
  int steps;
  int count;
  bool f_aum;
  bool f_maxmut;
  GraphI vac(jacta);
  GraphI new_one(jacta);
  
  
  for(int k = 0; k < num_exp; k++){
    f_maxmut = false;
    steps = 0;
    //Obtiene el genotipo fundador (red)
    ifstream inventario;
    bas.open_ifstream(inventario, "muestrario/" + bas.inttostring(seed) + "_inventario_2m.txt");
    ifstream muestra;
    bas.open_ifstream(muestra, "muestrario/" + bas.inttostring(seed) + "_muestra_2m.txt");
    
    red.get_dir_nw_from_file_bignw(nodes, muestra, inventario, k, samplesize);
    
    inventario.close();
    muestra.close();
    
    ocon << red.number_of_edges() << "\t";

    red.robustness_1mutation_cap(s, mgood, mbad, oc, ot, nodes, modules, predpar, rob);
    red.build_moma_d(matmod);
    orob << rob << "\t";
    omod << red.eval_mod(matmod, predpar) << "\t";
    
    vac.copy(red);
    while(!f_maxmut){
        f_aum = false;
        count = 0;
        vac.robustness_1mutation_cap(s, mgood, mbad, oc, ot, nodes, modules, predpar, last_rob);
        
        while(!f_aum && !f_maxmut){
            new_one.copy(vac);
            new_one.mutate(wol);
            new_one.interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
            if(fifu.cap_to_s(nodes, mgood, mbad, oc, ot, modules)){
                new_one.robustness_1mutation_cap(s, mgood, mbad, oc, ot, nodes, modules, predpar, new_rob);
                if(new_rob > last_rob){
                    vac.clear();
                    vac.copy(new_one);
                    f_aum = true;
                }
                else{
                    f_aum = false;
                    count++;
                }
            }
            new_one.clear();
            if(count == max_mut){
                f_maxmut = true;
            }
        }
        
        if(f_maxmut){
            vac.robustness_1mutation_cap(s, mgood, mbad, oc, ot, nodes, modules, predpar, rob);
            vac.build_moma_d(matmod);
            orob << rob << endl;;
            omod << vac.eval_mod(matmod, predpar) << endl;
            otime << steps << endl;
            break;
        }
        steps++;
    }
    
    ocon << vac.number_of_edges() << endl;
    
    vac.clear();
    red.clear();
  }
    
  
  ocon.close();
  orob.close();
  omod.close();
  otime.close();    
  
  return 0;

}
