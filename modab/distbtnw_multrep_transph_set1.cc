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

//funcion para conseguir distancia de saturacion
void distancia_max(Alea &jacta, int timetocheck, int crittosat, int sample_nt, int nugaps, int nodes, int *s, int *mgood, int *mbad, int **oc, int **ot, int modules, set<set<int>> &predpar, int ***intmods, double *tomod, double *frommod, double wol, int veces_pasos, int *dists){
  
  Basics bas(jacta);
  GraphI red(jacta); 
  GraphI red2(jacta);
  GraphI redtemp(jacta);
  FitnessI lafi(jacta);
  GraphI redtodist;
  bool f_dist, f_mod, f_cumple;
  double sum_dist = 0;
  int dist_sat, distcheck1, distcheck2;  
  int paso;
  int cuint;
  double modtemp;
  ifstream muestra;
  
  double **matmod;
  bas.create_array(matmod, nodes, nodes);
  for(int k = 0; k < sample_nt; k++){ 
      cuint = int(double(k)/100.0);
      f_dist = false;
      paso = 1;
      bas.open_ifstream(muestra,  "muestrario/" + bas.inttostring(k) + "/" + "nw_set1.txt");
      red.get_dir_nw_from_file(nodes, muestra);
        
      red.make_copy(red2);
      red.make_copy(redtodist);
      while(!f_dist){
          do{
              red2.make_copy(redtemp);
              redtemp.mutate_forced(wol);
              redtemp.interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
              f_cumple = lafi.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, nugaps);
              redtemp.build_moma_d(matmod);
              modtemp = redtemp.eval_mod(matmod, predpar);
              if(modtemp  > frommod[cuint] && modtemp < tomod[cuint])
                  f_mod = true;
              else
                  f_mod = false;
              if(f_cumple && f_mod){
                  red2.clear();
                  redtemp.make_copy(red2);
              }
              redtemp.clear();
          } while(!f_cumple || !f_mod);
          paso++;
          if(paso%timetocheck == 0){
              distcheck1 = red.distance_between_nw(red2);
              distcheck2 = red.distance_between_nw(redtodist);
              if(abs(distcheck1 - distcheck2) < crittosat)
                  f_dist = true;
              else{
                  redtodist.clear();
                  red2.make_copy(redtodist);
              }
          }
      }
      sum_dist += red.distance_between_nw(red2);
      red2.clear();
      red.clear();
      muestra.close();
      redtodist.clear();
    }
  
    dist_sat = (sum_dist/sample_nt);
    
    for(int i = 0; i < veces_pasos; i++)
        dists[i] = (dist_sat/veces_pasos)*(i+1);
    cout << "ter_dist" << endl;
    
    for(int n = 0; n < nodes; n++)
        delete matmod[n];
    delete[] matmod;
}

int main()
{
  
  int seed = 24125;
  int nodes = 24;
  int modules = 4;
  int sample_nt = 200;
  int inter_mod = 2;
  int at_step = 300;
  double prop_max_ac = 0.25;
  double wol = (double)3/(double)nodes;
  int genpermod = nodes/modules;
  int nugaps = modules;
  int rep_nw = 5;
  
  //parámetros de funcion distancia_max
  int timetocheck = 100;
  int crittosat = 20;
  
  double initial_mb = 0.1002;
  double initial_ma = 0.1934;
  double mues_sd = 0.0383591305190811;
  
  double *to, *from;
  to = new double[inter_mod];
  from = new double[inter_mod];

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI red(jacta);
  GraphI red2(jacta);
  GraphI redtemp(jacta);
  Cap_evolI cap(jacta);
  FitnessI lafi(jacta);

  int veces_pasos = 5;
  
  int veces_k = 1;
  double *k_mult;
  k_mult = new double[veces_k];
  string *kstr;
  kstr = new string[veces_k];
  
  for(int i = 0; i < veces_k; i++){
      k_mult[i] = double(i + 1);
      kstr[i] = bas.inttostring(i + 1) + "k";
  }
  
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
  
  int **ats;
  bas.create_array(ats, nugaps, nodes);
  cap.start_gaps_ics(ats, nugaps, modules, nodes, modgaps);
  
  int **mgood, **mbad;
  int ***oc, ***ot;
  
  bas.create_array(mgood, 2, nodes);
  bas.create_array(mbad, 2, nodes);
  bas.create_array(oc, 2, nodes, modules);
  bas.create_array(ot, 2, nodes, modules);
  
  int *dists;
  dists = new int[veces_pasos];
  
  int cuint;
  bool f_cumple, f_mod;
  
  int mean_count_pasos, count_pasos;
  double mean_tot_phen12, mean_tot_acphen12, mean_tot_phen1, mean_tot_acphen1, mean_tot_phen2, mean_tot_acphen2;
  int **tot_phen, **tot_acphen;
  bas.create_array(tot_phen, 3, nugaps);
  bas.create_array(tot_acphen, 3, nugaps);
  
  double mean_frac_newpht, mean_frac_newacph;
  double *frac_newpht, *frac_newacph;
  frac_newpht = new double[nugaps];
  frac_newacph = new double[nugaps];
  
  double mean_tot_acphens, mean_acphens_d1, mean_acphens_m1, mean_acphens_multchan, mean_sim_acphens, mean_mods_acphens, mean_frac_acphens_d1, mean_frac_acphens_m1;
  
  int *tot_acphens, *acphens_d1, *acphens_m1, *acphens_multchan;
  double *sim_acphens, *mods_acphens,  *frac_acphens_d1, *frac_acphens_m1;
  
  tot_acphens = new int[nugaps];
  acphens_d1 = new int[nugaps];
  acphens_m1 = new int[nugaps];
  acphens_multchan = new int[nugaps];
  sim_acphens = new double[nugaps];
  mods_acphens = new double[nugaps];
  frac_acphens_d1 = new double[nugaps];
  frac_acphens_m1 = new double[nugaps];
    
  double mean_tot_newpht, mean_newpht_g1, mean_newpht_d1, mean_newpht_m1, mean_newpht_multchan, mean_simnewpht, mean_modsnewpht, mean_frac_newpht_g1, mean_frac_newpht_d1, mean_frac_newpht_m1;
  
  double mean_tot_newacph, mean_newacph_g1, mean_newacph_d1, mean_newacph_m1, mean_newacph_multchan, mean_simnewacph, mean_modsnewacph, mean_frac_newacph_g1, mean_frac_newacph_d1, mean_frac_newacph_m1;
  
  int *tot_newpht, *newpht_g1, *newpht_d1, *newpht_m1, *newpht_multchan;
  double *simnewpht, *modsnewpht,  *frac_newpht_g1,  *frac_newpht_d1, *frac_newpht_m1;
  int *tot_newacph, *newacph_g1, *newacph_d1, *newacph_m1, *newacph_multchan;
  double *simnewacph, *modsnewacph, *frac_newacph_g1, *frac_newacph_d1, *frac_newacph_m1;
  
  tot_newpht = new int[nugaps];
  tot_newacph = new int[nugaps];
  newpht_g1 = new int[nugaps];
  newacph_g1 = new int[nugaps];
  newpht_d1 = new int[nugaps];
  newacph_d1 = new int[nugaps];
  newpht_m1 = new int[nugaps];
  newacph_m1 = new int[nugaps];
  newpht_multchan = new int[nugaps];
  newacph_multchan = new int[nugaps];
  
  simnewpht = new double[nugaps];
  simnewacph = new double[nugaps];
  modsnewpht = new double[nugaps];
  modsnewacph = new double[nugaps];
  frac_newpht_g1 = new double[nugaps];
  frac_newacph_g1 = new double[nugaps];
  frac_newpht_d1 = new double[nugaps];
  frac_newacph_d1 = new double[nugaps];
  frac_newpht_m1 = new double[nugaps];
  frac_newacph_m1 = new double[nugaps];
  
  
  int trunc;
  
  double **matmod;
  bas.create_array(matmod, nodes, nodes);
  double mean_mod1, mean_mod2, mod1, mod2, modtemp, mean_edges_red2, edges_red2;
  
  double dist_btw_nw; 
  
  ofstream **salida;
  salida = new ofstream *[inter_mod];
  for (int i = 0; i < inter_mod; i++)
      salida[i] = new ofstream[3];
  
  ifstream muestra;
  string *mods;
  mods = new string[inter_mod];
  
  for(int l = 0; l < veces_k; l++){
  
    mods[0] = "modb_" + kstr[l];
    mods[1] = "moda_" + kstr[l];
    
    to[0] = initial_mb + k_mult[l]*mues_sd;
    from[0] = initial_mb - k_mult[l]*mues_sd;
    to[1] = initial_ma + k_mult[l]*mues_sd;
    from[1] = initial_ma - k_mult[l]*mues_sd;
    
    distancia_max(jacta, timetocheck, crittosat, sample_nt, nugaps, nodes, s, mgood[0], mbad[0], oc[0], ot[0], modules, predpar, intmods, to, from, wol, veces_pasos, dists);
    
    for(int i = 0; i < inter_mod; i++){
        bas.open_ofstream(salida[i][0],  "distbtwnw/" + bas.inttostring(seed) + "_" + mods[i] + "_" + bas.inttostring(rep_nw) + "rep_general_set1.txt");
        salida[i][0] << "nu_pasos" << "\t" << "tot_phen" << "\t" << "tot_phen1" << "\t" << "tot_phen2" << "\t" << "tot_acphen" << "\t" << "tot_acphen1" << "\t" << "tot_acphen2" << "\t" << "frac_newpht" << "\t" << "frac_newacph"  << "\t" << "mod_red1" << "\t" << "mod_red2" << "\t" << "edges_red2" << "\t" << "dist_btw_nw" << "\t" << "Trunc" << endl;
        
        bas.open_ofstream(salida[i][1],  "distbtwnw/" + bas.inttostring(seed) + "_" + mods[i] + "_" + bas.inttostring(rep_nw) + "rep_acphens_set1.txt");
        salida[i][1] << "tot_acphens" << "\t" << "sim_acphens" << "\t" << "mods_acphens" << "\t" << "acphens_d1" << "\t" << "acphens_m1" << "\t" << "frac_acphens_d1" << "\t" << "frac_acphens_m1" << "\t" << "acphens_multchan" << "\t" << "dist_btw_nw"  << endl;
        
        bas.open_ofstream(salida[i][2],  "distbtwnw/" + bas.inttostring(seed) + "_" + mods[i] + "_" + bas.inttostring(rep_nw) + "rep_newphen_set1.txt");
        salida[i][2] << "tot_newpht" << "\t" << "simnewpht" << "\t" << "modsnewpht" << "\t" << "newpht_g1" << "\t" << "newpht_d1" << "\t" << "newpht_m1" << "\t" << "frac_newpht_g1" << "\t" << "frac_newpht_d1" << "\t" << "frac_newpht_m1" << "\t" << "newpht_multchan" << "\t" << "tot_newacph" << "\t" << "simnewacph" << "\t" << "modsnewacph" << "\t" << "newacph_g1" << "\t" << "newacph_d1" << "\t" << "newacph_m1" << "\t" << "frac_newacph_g1" << "\t" << "frac_newacph_d1" << "\t" << "frac_newacph_m1" << "\t" << "newacph_multchan"  << "\t" << "dist_btw_nw" << endl;
    
    }
    for(int i = 0; i < sample_nt; i++){
        cuint = int(double(i)/100.0);
        cout << cuint <<  " " << i << endl;
        bas.open_ifstream(muestra,  "muestrario/" + bas.inttostring(i) + "/" + "nw_set1.txt");
        
        red.get_dir_nw_from_file(nodes, muestra);
        
        red.build_moma_d(matmod);
        mod1 = red.eval_mod(matmod, predpar);
        if(!(mod1  > from[cuint] && mod1 < to[cuint])){
            cout << "[Error]: Graph " << i << " do not have desired modularity.\n";
            exit(1);
        }
    //         cout << mod1 << endl;
        for(int j = 0; j < veces_pasos; j++){
            mean_count_pasos = 0;
            mean_tot_phen12 = 0;
            mean_tot_acphen12 = 0;
            mean_tot_phen1 = 0;
            mean_tot_acphen1 = 0;
            mean_tot_phen2 = 0;
            mean_tot_acphen2 = 0;
            mean_mod1 = 0;
            mean_mod2 = 0;
            mean_edges_red2 = 0;
            
            mean_frac_newpht = 0;
            mean_frac_newacph = 0;
            
            mean_tot_acphens = 0;
            mean_sim_acphens = 0;
            mean_mods_acphens = 0;
            mean_acphens_d1 = 0;
            mean_acphens_m1 = 0;
            mean_frac_acphens_d1 = 0;
            mean_frac_acphens_m1 = 0;
            mean_acphens_multchan = 0;
            
            mean_tot_newpht = 0;
            mean_simnewpht = 0;
            mean_modsnewpht = 0;
            mean_newpht_g1 = 0;
            mean_newpht_d1 = 0;
            mean_newpht_m1 = 0;
            mean_frac_newpht_g1 = 0;
            mean_frac_newpht_d1 = 0;
            mean_frac_newpht_m1 = 0;
            mean_newpht_multchan = 0;
            
            mean_tot_newacph = 0;
            mean_simnewacph = 0;
            mean_modsnewacph = 0;
            mean_newacph_g1 = 0;
            mean_newacph_d1 = 0;
            mean_newacph_m1 = 0;
            mean_frac_newacph_g1 = 0;
            mean_frac_newacph_d1 = 0;
            mean_frac_newacph_m1 = 0;
            mean_newacph_multchan = 0;
            
            for(int nurep = 0; nurep < rep_nw; nurep++){
                bas.fillmat0(tot_phen, 3, nugaps);
                bas.fillmat0(tot_acphen, 3, nugaps);
                
                bas.fillv0(frac_newpht, nugaps);
                bas.fillv0(frac_newacph, nugaps);
                bas.fillv0(tot_acphens, nugaps);
                bas.fillv0(sim_acphens, nugaps);
                bas.fillv0(mods_acphens, nugaps);
                bas.fillv0(acphens_d1, nugaps);
                bas.fillv0(acphens_m1, nugaps);
                bas.fillv0(acphens_multchan, nugaps);
                bas.fillv0(frac_acphens_d1, nugaps);
                bas.fillv0(frac_acphens_m1, nugaps);
                
                bas.fillv0(tot_newpht, nugaps);
                bas.fillv0(tot_newacph, nugaps);
                bas.fillv0(newpht_g1, nugaps);
                bas.fillv0(newpht_d1, nugaps);
                bas.fillv0(simnewpht, nugaps);
                bas.fillv0(simnewacph, nugaps);
                bas.fillv0(modsnewpht, nugaps);
                bas.fillv0(modsnewacph, nugaps);
                bas.fillv0(newacph_g1, nugaps);
                bas.fillv0(newacph_d1, nugaps);
                bas.fillv0(newpht_m1, nugaps);
                bas.fillv0(newacph_m1, nugaps);
                bas.fillv0(newpht_multchan, nugaps);
                bas.fillv0(newacph_multchan, nugaps);
                
                bas.fillv0(frac_newpht_g1, nugaps);
                bas.fillv0(frac_newacph_g1, nugaps);
                bas.fillv0(frac_newpht_d1, nugaps);
                bas.fillv0(frac_newacph_d1, nugaps);
                bas.fillv0(frac_newpht_m1, nugaps);
                bas.fillv0(frac_newacph_m1, nugaps);
                
                red.make_copy(red2);
                count_pasos = 0;
                do{
                    do{
                        red2.make_copy(redtemp);
                        redtemp.mutate_forced(wol);
                        redtemp.interactions_to_attractor(s, mgood[1], mbad[1], oc[1], ot[1], nodes, modules, predpar);
                        f_cumple = lafi.cap_to_s_gaps(nodes, mgood[1], mbad[1], oc[1], ot[1], modules, intmods, nugaps);
                        redtemp.build_moma_d(matmod);
                        modtemp = redtemp.eval_mod(matmod, predpar);
                        if(modtemp  > from[cuint] && modtemp < to[cuint])
                            f_mod = true;
                        else
                            f_mod = false;
                        if(f_cumple && f_mod){
                            red2.clear();
                            redtemp.make_copy(red2);
                        }
                        redtemp.clear();
                    } while(!f_cumple || !f_mod);
                    count_pasos++;
                    dist_btw_nw = red.distance_between_nw(red2);
                }while(dist_btw_nw < dists[j]);
                
                red.access_phen_translape_gaps(red2, ats, nugaps, s, mgood, mbad, oc, ot, nodes, modules, predpar, intmods, at_step, prop_max_ac, tot_phen, tot_acphen, tot_acphens, sim_acphens, mods_acphens, acphens_d1, acphens_m1, acphens_multchan, tot_newpht, simnewpht, modsnewpht, newpht_g1, newpht_d1, newpht_m1, newpht_multchan, tot_newacph, simnewacph, modsnewacph, newacph_g1, newacph_d1, newacph_m1, newacph_multchan, trunc);
                
                dist_btw_nw = red.distance_between_nw(red2);
                
                red2.build_moma_d(matmod);
                mod2 = red2.eval_mod(matmod, predpar);
                edges_red2 = red2.number_of_edges();
                red2.clear();
            
                for(int k = 0; k < nugaps; k++){
                    tot_phen[2][k] = tot_phen[0][k] + tot_phen[1][k];
                    if(tot_phen[1][k] != 0)
                        frac_newpht[k] = double(tot_newpht[k])/double(tot_phen[1][k]);
                    
                    tot_acphen[2][k] = tot_acphen[0][k] + tot_acphen[1][k];
                    if(tot_acphen[1][k] != 0)
                        frac_newacph[k] = double(tot_newacph[k])/double(tot_acphen[1][k]);
                }
                
                mean_count_pasos += count_pasos;
                mean_tot_phen12 += bas.get_mean(tot_phen[2], nugaps);
                mean_tot_phen1 += bas.get_mean(tot_phen[0], nugaps);
                mean_tot_phen2 += bas.get_mean(tot_phen[1], nugaps);
                mean_tot_acphen12 += bas.get_mean(tot_phen[2], nugaps);
                mean_tot_acphen1 += bas.get_mean(tot_phen[0], nugaps);
                mean_tot_acphen2 += bas.get_mean(tot_phen[1], nugaps);
                mean_frac_newpht += bas.get_mean(frac_newpht, nugaps);
                mean_frac_newacph += bas.get_mean(frac_newacph, nugaps);
                mean_mod1 += mod1;
                mean_mod2 += mod2;
                mean_edges_red2 += edges_red2;
                
                for(int k = 0; k < nugaps; k++){
                    if(tot_acphens[k] != 0)
                        frac_acphens_d1[k] = double(acphens_d1[k])/double(tot_acphens[k]);
                    
                    if(acphens_multchan[k] != 0)
                        frac_acphens_m1[k] = double(acphens_m1[k])/double(acphens_multchan[k]);
                }
                
                mean_tot_acphens += bas.get_mean(tot_acphens, nugaps);
                mean_sim_acphens += bas.get_mean(sim_acphens, nugaps);
                mean_mods_acphens += bas.get_mean(mods_acphens, nugaps);
                mean_acphens_d1 += bas.get_mean(acphens_d1, nugaps);
                mean_acphens_m1 += bas.get_mean(acphens_m1, nugaps);
                mean_frac_acphens_d1 += bas.get_mean(frac_acphens_d1, nugaps);
                mean_frac_acphens_m1 += bas.get_mean(frac_acphens_m1, nugaps);
                mean_acphens_multchan += bas.get_mean(acphens_multchan, nugaps);
                
                for(int k = 0; k < nugaps; k++){
                    if(tot_newacph[k] != 0){
                        frac_newpht_g1[k] = double(newpht_g1[k])/double(tot_newpht[k]);
                        frac_newacph_g1[k] = double(newacph_g1[k])/double(tot_newacph[k]);
                        
                        frac_newpht_d1[k] = double(newpht_d1[k])/double(tot_newpht[k]);
                        frac_newacph_d1[k] = double(newacph_d1[k])/double(tot_newacph[k]);
                    }
                    
                    if(newacph_multchan[k] != 0){
                        frac_newpht_m1[k] = double(newpht_m1[k])/double(newpht_multchan[k]);
                        frac_newacph_m1[k] = double(newacph_m1[k])/double(newacph_multchan[k]);
                    }
                }
                mean_tot_newpht += bas.get_mean(tot_newpht, nugaps);
                mean_simnewpht += bas.get_mean(simnewpht, nugaps);
                mean_modsnewpht += bas.get_mean(modsnewpht, nugaps);
                mean_newpht_g1 += bas.get_mean(newpht_g1, nugaps);
                mean_newpht_d1 += bas.get_mean(newpht_d1, nugaps);
                mean_newpht_m1 += bas.get_mean(newpht_m1, nugaps);
                mean_frac_newpht_g1 += bas.get_mean(frac_newpht_g1, nugaps);
                mean_frac_newpht_d1 += bas.get_mean(frac_newpht_d1, nugaps);
                mean_frac_newpht_m1 += bas.get_mean(frac_newpht_m1, nugaps);
                mean_newpht_multchan += bas.get_mean(newpht_multchan, nugaps);
                
                mean_tot_newacph += bas.get_mean(tot_newacph, nugaps);
                mean_simnewacph += bas.get_mean(simnewacph, nugaps);
                mean_modsnewacph += bas.get_mean(modsnewacph, nugaps);
                mean_newacph_g1 += bas.get_mean(newacph_g1, nugaps);
                mean_newacph_d1 += bas.get_mean(newacph_d1, nugaps);
                mean_newacph_m1 += bas.get_mean(newacph_m1, nugaps);
                mean_frac_newacph_g1 += bas.get_mean(frac_newacph_g1, nugaps);
                mean_frac_newacph_d1 += bas.get_mean(frac_newacph_d1, nugaps);
                mean_frac_newacph_m1 += bas.get_mean(frac_newacph_m1, nugaps);
                mean_newacph_multchan += bas.get_mean(newacph_multchan, nugaps);
                
            }
            
            salida[cuint][0] << mean_count_pasos/double(rep_nw) << "\t" << mean_tot_phen12/double(rep_nw) << "\t" << mean_tot_phen1/double(rep_nw) << "\t" << mean_tot_phen2/double(rep_nw) << "\t" << mean_tot_acphen12/double(rep_nw) << "\t" << mean_tot_acphen1/double(rep_nw) << "\t" << mean_tot_acphen2/double(rep_nw) << "\t" << mean_frac_newpht/double(rep_nw) << "\t" << mean_frac_newacph/double(rep_nw) << "\t" << mean_mod1/double(rep_nw) << "\t" << mean_mod2/double(rep_nw) << "\t" << mean_edges_red2/double(rep_nw) << "\t" << dist_btw_nw << "\t" << trunc << endl;
            
            salida[cuint][1] << mean_tot_acphens/double(rep_nw) << "\t" << mean_sim_acphens/double(rep_nw) << "\t" << mean_mods_acphens/double(rep_nw) << "\t" << mean_acphens_d1/double(rep_nw) << "\t" << mean_acphens_m1/double(rep_nw) << "\t" << mean_frac_acphens_d1/double(rep_nw) << "\t" << mean_frac_acphens_m1/double(rep_nw) << "\t" << mean_acphens_multchan/double(rep_nw) << "\t" << dist_btw_nw << endl;
            
            salida[cuint][2] << mean_tot_newpht/double(rep_nw) << "\t" << mean_simnewpht/double(rep_nw) << "\t" << mean_modsnewpht/double(rep_nw) << "\t" << mean_newpht_g1/double(rep_nw) << "\t" << mean_newpht_d1/double(rep_nw) << "\t" << mean_newpht_m1/double(rep_nw) << "\t" << mean_frac_newpht_g1/double(rep_nw) << "\t" << mean_frac_newpht_d1/double(rep_nw) << "\t" << mean_frac_newpht_m1/double(rep_nw) << "\t" << mean_newpht_multchan/double(rep_nw) << "\t";
            
            salida[cuint][2] << mean_tot_newacph/double(rep_nw) << "\t" << mean_simnewacph/double(rep_nw) << "\t" << mean_modsnewacph/double(rep_nw) << "\t" << mean_newacph_g1/double(rep_nw) << "\t" << mean_newacph_d1/double(rep_nw) << "\t" << mean_newacph_m1/double(rep_nw) << "\t" << mean_frac_newacph_g1/double(rep_nw) << "\t" << mean_frac_newacph_d1/double(rep_nw) << "\t" << mean_frac_newacph_m1/double(rep_nw)<< "\t" << mean_newacph_multchan/double(rep_nw) << "\t" << dist_btw_nw << endl;
        }
        
        red.clear();
        
        muestra.close();
    }

    for(int i = 0; i < inter_mod; i++){
        salida[i][0].close();
        salida[i][1].close();
        salida[i][2].close();
    }
    
   
  }
  
  delete[] s;
  jacta.close_rng();
  cout << "ya" << endl;
  return 0;

}
