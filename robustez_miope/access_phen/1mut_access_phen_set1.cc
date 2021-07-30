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
  int nugaps = modules;
  int genpermod = nodes/modules;   
  double prop_max_ac = 0.5;

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI red(jacta);
  Cap_evolI chw(jacta);
  FitnessI lafi(jacta);

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

  
  int *s;
  s = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       if ((i%2)==0)
            s[i] = 1;
       else
            s[i] = -1;
  }
  
  int **ats;
  bas.create_array(ats, nugaps, nodes);
  chw.start_gaps_ics(ats, nugaps, modules, nodes, modgaps);
  
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
  
  int *tot_pht, *pht_g1, *pht_d1, *pht_10, *pht_m1, *pht_multchan;
  int *tot_acph, *acph_g1, *acph_d1, *acph_10, *acph_m1, *acph_multchan;
  double *simpht, *modspht, *frac_pht_g1, *frac_pht_d1, *frac_pht_10, *frac_pht_m1;
  double *simacph, *modsacph, *frac_acph_g1, *frac_acph_d1, *frac_acph_10, *frac_acph_m1;
  int nu_multchan;
  
  tot_pht = new int[nugaps];
  tot_acph = new int[nugaps];
  
  pht_d1 = new int[nugaps];
  acph_d1 = new int[nugaps];
  
  pht_g1 = new int[nugaps];
  acph_g1 = new int[nugaps];
  
  pht_10 = new int[nugaps];
  acph_10 = new int[nugaps];
  
  pht_m1 = new int[nugaps];
  acph_m1 = new int[nugaps];
  
  pht_multchan = new int[nugaps];
  acph_multchan = new int[nugaps];
  
  simpht = new double[nugaps];
  simacph = new double[nugaps];
  
  modspht = new double[nugaps];
  modsacph = new double[nugaps];
  
  frac_pht_d1 = new double[nugaps];
  frac_acph_d1 = new double[nugaps];
  
  frac_pht_g1 = new double[nugaps];
  frac_acph_g1 = new double[nugaps];
  
  frac_pht_10 = new double[nugaps];
  frac_acph_10 = new double[nugaps];

  frac_pht_m1 = new double[nugaps];
  frac_acph_m1 = new double[nugaps];
  
  
  double dist_acph, phen_size;
  int maxsize, trunc;
  
  string directorio = "muestrario";
  
  ifstream muestra;
  ifstream inventario;
  ofstream orob;
  
  bas.open_ofstream(orob, "access_phen/"  + bas.inttostring(seed) + "_allmut_set1.txt");
  orob << "tot_pht" << "\t" << "simpht" << "\t" << "modspht" << "\t" << "pht_g1" << "\t" << "pht_d1" << "\t" << "pht_10" << "\t" << "pht_m1" << "\t" << "frac_pht_g1" << "\t" << "frac_pht_d1" << "\t" << "frac_pht_10" << "\t" << "frac_pht_m1" << "\t" << "tot_acph" << "\t" << "simacph" << "\t" << "modsacph" << "\t" << "acph_g1" << "\t" << "acph_d1" << "\t" << "acph_10" << "\t" << "acph_m1" << "\t" << "frac_acph_g1" << "\t" << "frac_acph_d1" << "\t" << "frac_acph_10" << "\t" << "frac_acph_m1" << "\t" << "att_max_size" << "\t" << "att_mean_size"  << "\t" << "dist_btw_acph" << "\t" << "Trunc" << endl;
  //Obtiene los genotipos extremos
  for (int i = 0; i < nu_exp; i++){
//       cout << i << endl;
      bas.open_ifstream(muestra, directorio + "/" + bas.inttostring(seed) + "_muestra_set1.txt");
      bas.open_ifstream(inventario, directorio + "/" + bas.inttostring(seed) + "_inventario_set1.txt");

      red.get_dir_nw_from_file_bignw(nodes, muestra, inventario, i, sample_nt);
      
      bas.fillv0(tot_pht, nugaps);
      bas.fillv0(tot_acph, nugaps);
      bas.fillv0(pht_g1, nugaps);
      bas.fillv0(acph_g1, nugaps);
      bas.fillv0(pht_d1, nugaps);
      bas.fillv0(acph_d1, nugaps);
      bas.fillv0(pht_10, nugaps);
      bas.fillv0(acph_10, nugaps);
      bas.fillv0(pht_m1, nugaps);
      bas.fillv0(acph_m1, nugaps);
      bas.fillv0(pht_multchan, nugaps);
      bas.fillv0(acph_multchan, nugaps);
      bas.fillv0(simpht, nugaps);
      bas.fillv0(simacph, nugaps);
      bas.fillv0(modspht, nugaps);
      bas.fillv0(modsacph, nugaps);
      bas.fillv0(frac_pht_g1, nugaps);
      bas.fillv0(frac_acph_g1, nugaps);
      bas.fillv0(frac_pht_d1, nugaps);
      bas.fillv0(frac_acph_d1, nugaps);
      bas.fillv0(frac_pht_10, nugaps);
      bas.fillv0(frac_acph_10, nugaps);
      bas.fillv0(frac_pht_m1, nugaps);
      bas.fillv0(frac_acph_m1, nugaps);
      nu_multchan = 0;
      
      red.access_phen_1mutation_gaps(ats, nugaps, s, mgood, mbad, oc, ot, nodes, modules, predpar, intmods, at_step, prop_max_ac, tot_pht, simpht, modspht, pht_g1, pht_d1, pht_10, pht_m1, pht_multchan, tot_acph, simacph, modsacph, acph_g1, acph_d1, acph_10, acph_m1, acph_multchan, maxsize, phen_size, dist_acph, trunc);
      for(int ic = 0; ic < nugaps; ic++){
          if(tot_acph[ic] != 0){
            frac_pht_g1[ic] = pht_g1[ic]/double(tot_pht[ic]);
            frac_pht_d1[ic] = pht_d1[ic]/double(tot_pht[ic]);
            frac_pht_10[ic] = pht_10[ic]/double(tot_pht[ic]);
            
            frac_acph_g1[ic] = acph_g1[ic]/double(tot_acph[ic]);
            frac_acph_d1[ic] = acph_d1[ic]/double(tot_acph[ic]);
            frac_acph_10[ic] = acph_10[ic]/double(tot_acph[ic]);
          }

          if(acph_multchan[ic] != 0){
            frac_pht_m1[ic] = pht_m1[ic]/double(pht_multchan[ic]);
            frac_acph_m1[ic] = acph_m1[ic]/double(acph_multchan[ic]);
            nu_multchan++;
          }
      }
    
      orob << bas.get_mean(tot_pht, nugaps) << "\t" << bas.get_mean(simpht, nugaps) << "\t" << double(bas.sumatoria(modspht, nugaps))/double(nu_multchan) << "\t" << bas.get_mean(pht_g1, nugaps) << "\t" << bas.get_mean(pht_d1, nugaps) << "\t" << bas.get_mean(pht_10, nugaps)  << "\t" << double(bas.sumatoria(pht_m1, nugaps))/double(nu_multchan) << "\t" << bas.get_mean(frac_pht_g1, nugaps) << "\t" << bas.get_mean(frac_pht_d1, nugaps) << "\t" << bas.get_mean(frac_pht_10, nugaps) << "\t" << double(bas.sumatoria(frac_pht_m1, nugaps))/double(nu_multchan) << "\t"; 
    
      orob << bas.get_mean(tot_acph, nugaps) << "\t" << bas.get_mean(simacph, nugaps) << "\t" << double(bas.sumatoria(modsacph, nugaps))/double(nu_multchan) << "\t" << bas.get_mean(acph_g1, nugaps) << "\t" << bas.get_mean(acph_d1, nugaps) << "\t" << bas.get_mean(acph_10, nugaps) << "\t" << double(bas.sumatoria(acph_m1, nugaps))/double(nu_multchan) << "\t" << bas.get_mean(frac_acph_g1, nugaps) << "\t" << bas.get_mean(frac_acph_d1, nugaps)  << "\t" << bas.get_mean(frac_acph_10, nugaps) << "\t" << double(bas.sumatoria(frac_acph_m1, nugaps))/double(nu_multchan) << "\t";
    
      orob << maxsize << "\t" << phen_size << "\t" << dist_acph << "\t" << trunc << endl; 
      
      red.clear();
      
      muestra.close();
      inventario.close();
  }
  orob.close();
  
  jacta.close_rng();
  cout << "ya" << endl;  
}
