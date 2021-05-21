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
  int genpermod = nodes/modules;
  int nugaps = modules;
  
  int max_mut = 2*nodes;
  double wol = double(conectividad)/double(nodes*nodes);
     
  Alea jacta(seed);
  Basics bas(jacta);
  Cap_evolI cap(jacta);
  GraphI red(jacta);
  FitnessI fifu(jacta);
  
  int** modgaps;
  bas.create_array(modgaps, nugaps, modules);
  for(int j = 0; j < nugaps; j++)
      for(int k = 0; k < modules; k++){
          if(j == k)
              modgaps[j][k] = 1;
          else
              modgaps[j][k] = 0;
      }

      
  int ***intmods;
  bas.create_array(intmods, nodes, nugaps, modules);
  
  for(int i = 0; i < nodes; i++){
      for(int j = 0; j < nugaps; j++){
          for(int k = 0; k < modules; k++){
              if ((int)(i/genpermod) == k)
                  intmods[i][j][k] = 0;
              else{
                  if (j == 0)
                      intmods[i][j][k] = 1;
                  else
                      if((modgaps[0][(int)(i/genpermod)] == modgaps[j][(int)(i/genpermod)]) == (modgaps[0][k] == modgaps[j][k])) 
                          intmods[i][j][k] = 1;
                      else 
                          intmods[i][j][k] = -1;
              }
          }
      }
  }
  
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
  double last_mod, new_mod;
  
  ofstream ocon;
  bas.open_ofstream(ocon, "miope_mod/" +  bas.inttostring(seed) + "_con_set1.txt");
  ocon << "Initial" << "\t" << "Final" << endl;
    
  ofstream omod;
  bas.open_ofstream(omod, "miope_mod/" + bas.inttostring(seed) + "_mod_set1.txt");
  omod << "Initial" << "\t" << "Final" << endl;
  
  ofstream orob;
  bas.open_ofstream(orob, "miope_mod/" + bas.inttostring(seed) + "_rob_set1.txt");
  orob << "Initial" << "\t" << "Final" << endl;
  
  ofstream otime;
  bas.open_ofstream(otime, "miope_mod/" + bas.inttostring(seed) + "_time_set1.txt");
  
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
    bas.open_ifstream(inventario, "muestrario/" + bas.inttostring(seed) + "_inventario_set1.txt");
    ifstream muestra;
    bas.open_ifstream(muestra, "muestrario/" + bas.inttostring(seed) + "_muestra_set1.txt");
    
    red.get_dir_nw_from_file_bignw(nodes, muestra, inventario, k, samplesize);
    
    inventario.close();
    muestra.close();
    
    ocon << red.number_of_edges() << "\t";

    if(!red.robustness_1mutation_gaps(s, mgood, mbad, oc, ot, nodes, modules, predpar, intmods, nugaps, rob)){
        cout << "[Error]: Network did not produce the modular attractors.\n";
        exit(1);
    }
    red.build_moma_d(matmod);
    orob << rob << "\t";
    omod << red.eval_mod(matmod, predpar) << "\t";
    
    vac.copy(red);
    while(!f_maxmut){
        f_aum = false;
        count = 0;
        vac.build_moma_d(matmod);
        last_mod = vac.eval_mod(matmod, predpar);
        
        while(!f_aum && !f_maxmut){
            new_one.copy(vac);
            new_one.mutate(wol);
            new_one.interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
            if(fifu.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, nugaps)){
                new_one.build_moma_d(matmod);
                new_mod = new_one.eval_mod(matmod, predpar);
                if(new_mod > last_mod){
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
            vac.robustness_1mutation_gaps(s, mgood, mbad, oc, ot, nodes, modules, predpar, intmods, nugaps, rob);
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
