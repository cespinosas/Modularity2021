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
  int nuic = pow(2, (modules-1));
  double prop_max_ac = 0.5;

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

  
  int *tot_pht, *pht_g1, *pht_d1, *pht_m1, *pht_multchan;
  int *tot_acph, *acph_g1, *acph_d1, *acph_m1, *acph_multchan;
  double *simpht, *modspht,  *frac_pht_g1,  *frac_pht_d1, *frac_pht_m1;
  double *simacph, *modsacph, *frac_acph_g1, *frac_acph_d1, *frac_acph_m1;
  int nu_multchan;
  
  tot_pht = new int[nuic];
  tot_acph = new int[nuic];
  
  pht_g1 = new int[nuic];
  acph_g1 = new int[nuic];
  
  pht_d1 = new int[nuic];
  acph_d1 = new int[nuic];
  
  pht_m1 = new int[nuic];
  acph_m1 = new int[nuic];
  
  pht_multchan = new int[nuic];
  acph_multchan = new int[nuic];
  
  simpht = new double[nuic];
  simacph = new double[nuic];
  
  modspht = new double[nuic];
  modsacph = new double[nuic];
  
  frac_pht_g1 = new double[nuic];
  frac_acph_g1 = new double[nuic];
  
  frac_pht_d1 = new double[nuic];
  frac_acph_d1 = new double[nuic];
  
  frac_pht_m1 = new double[nuic];
  frac_acph_m1 = new double[nuic];
  
  double dist_acph, phen_size;
  int maxsize, trunc;
  
  string directorio = "muestrario";
  
  ifstream muestra;
  ifstream inventario;
  ofstream orob;
  
  bas.open_ofstream(orob, "access_phen/"  + bas.inttostring(seed) + "_allmut_2m.txt");
  orob << "tot_pht" << "\t" << "simpht" << "\t" << "modspht" << "\t" << "pht_g1" << "\t" << "pht_d1" << "\t" << "pht_m1" << "\t" << "frac_pht_g1" << "\t" << "frac_pht_d1" << "\t" << "frac_pht_m1" << "\t" << "tot_acph" << "\t" << "simacph" << "\t" << "modsacph" << "\t" << "acph_g1" << "\t" << "acph_d1" << "\t" << "acph_m1" << "\t" << "frac_acph_g1" << "\t" << "frac_acph_d1" << "\t" << "frac_acph_m1" << "\t" << "att_max_size" << "\t" << "att_mean_size"  << "\t" << "dist_btw_acph" << "\t" << "Trunc" << endl;
  //Obtiene los genotipos extremos
  for (int i = 0; i < nu_exp; i++){
//       cout << i << endl;
      bas.open_ifstream(muestra, directorio + "/" + bas.inttostring(seed) + "_muestra_2m.txt");
      bas.open_ifstream(inventario, directorio + "/" + bas.inttostring(seed) + "_inventario_2m.txt");

      red.get_dir_nw_from_file_bignw(nodes, muestra, inventario, i, sample_nt);
      
      bas.fillv0(tot_pht, nuic);
      bas.fillv0(tot_acph, nuic);
      bas.fillv0(pht_g1, nuic);
      bas.fillv0(acph_g1, nuic);
      bas.fillv0(pht_d1, nuic);
      bas.fillv0(acph_d1, nuic);
      bas.fillv0(pht_m1, nuic);
      bas.fillv0(acph_m1, nuic);
      bas.fillv0(pht_multchan, nuic);
      bas.fillv0(acph_multchan, nuic);
      bas.fillv0(simpht, nuic);
      bas.fillv0(simacph, nuic);
      bas.fillv0(modspht, nuic);
      bas.fillv0(modsacph, nuic);
      bas.fillv0(frac_pht_g1, nuic);
      bas.fillv0(frac_acph_g1, nuic);
      bas.fillv0(frac_pht_d1, nuic);
      bas.fillv0(frac_acph_d1, nuic);
      bas.fillv0(frac_pht_m1, nuic);
      bas.fillv0(frac_acph_m1, nuic);
      nu_multchan = 0;
      
      red.access_phen_1mutation_cap(ats, nuic, s, mgood, mbad, oc, ot, nodes, modules, predpar, at_step, prop_max_ac, tot_pht, simpht, modspht, pht_g1, pht_d1, pht_m1, pht_multchan, tot_acph, simacph, modsacph, acph_g1, acph_d1, acph_m1, acph_multchan, maxsize, phen_size, dist_acph, trunc);
      for(int ic = 0; ic < nuic; ic++){

          if(tot_acph[ic] != 0){
            frac_pht_g1[ic] = pht_g1[ic]/double(tot_pht[ic]);
            frac_acph_g1[ic] = acph_g1[ic]/double(tot_acph[ic]);
            
            frac_pht_d1[ic] = pht_d1[ic]/double(tot_pht[ic]);
            frac_acph_d1[ic] = acph_d1[ic]/double(tot_acph[ic]);
          }
          
          if(acph_multchan[ic] != 0){
            frac_pht_m1[ic] = pht_m1[ic]/double(pht_multchan[ic]);
            frac_acph_m1[ic] = acph_m1[ic]/double(acph_multchan[ic]);
            nu_multchan++;
          }
      }
    
      orob << bas.get_mean(tot_pht, nuic) << "\t" << bas.get_mean(simpht, nuic) << "\t" << double(bas.sumatoria(modspht, nuic))/double(nu_multchan) << "\t" << bas.get_mean(pht_g1, nuic) << "\t" << bas.get_mean(pht_d1, nuic)  << "\t" << double(bas.sumatoria(pht_m1, nuic))/double(nu_multchan) << "\t" << bas.get_mean(frac_pht_g1, nuic) << "\t" << bas.get_mean(frac_pht_d1, nuic) << "\t" << double(bas.sumatoria(frac_pht_m1, nuic))/double(nu_multchan) << "\t"; 
    
      orob << bas.get_mean(tot_acph, nuic) << "\t" << bas.get_mean(simacph, nuic) << "\t" << double(bas.sumatoria(modsacph, nuic))/double(nu_multchan) << "\t" << bas.get_mean(acph_g1, nuic) << "\t" << bas.get_mean(acph_d1, nuic) << "\t" << double(bas.sumatoria(acph_m1, nuic))/double(nu_multchan) << "\t" << bas.get_mean(frac_acph_g1, nuic) << "\t" << bas.get_mean(frac_acph_d1, nuic) << "\t" << double(bas.sumatoria(frac_acph_m1, nuic))/double(nu_multchan) << "\t";
    
      orob << maxsize << "\t" << phen_size << "\t" << dist_acph << "\t" << trunc << endl; 
      
      red.clear();
      
      muestra.close();
      inventario.close();
  }
  orob.close();
  
  jacta.close_rng();
  cout << "ya" << endl;  
}

