#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include <string>
#include "basics.h"
#include "graphi.h"
#include "fitnessi.h"
#include "cap_evoli.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

Cap_evolI::Cap_evolI()
{
  yaparamsmu = false;
  yapart_mues = false;
}

Cap_evolI::Cap_evolI(Alea& jacta)
{
  start_rng(jacta);
  yaparamsmu = false;
  yapart_mues = false;
}

void Cap_evolI::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

void Cap_evolI::set_params_muestreo(int samplesizep, int nugenp, int stepsizep, int conectividadp, int modulesp)
{
  samplesize = samplesizep;
  nugen = nugenp;
  stepsize = stepsizep;
  conectividad = conectividadp;
  modules = modulesp;
  yaparamsmu = true;
}

void Cap_evolI::set_partition(set<set<int> > &predparp)
{
    if(modules != (int)predparp.size()){
      cout << "[Error]: Number of modules and modules in given partition are not the same in Cap_evolI::set_partition.\n";
      exit(1);
    }
  set<set<int> >::iterator it;
  for(it = predparp.begin(); it != predparp.end(); it++)
    predpar.insert(*it);
  
  modkey = new set<set<int> >::iterator [predpar.size()];
  
  for(int i = 0; i < modules; i++){
      modkey[i] = predpar.begin();
      for(int j = 0; j < i; j++)
          modkey[i]++;
  }
  yapart_mues = true;
  
}


void Cap_evolI::print_params_mues(ostream& sal)
{
  sal << "Integer_weights\n";
  sal << "Number_of_genes: " << nugen << endl;
  sal << "Sample_size: " << samplesize<< endl;
  sal << "Sample_connectivity: " << conectividad << endl;
  sal << "Stepsize: " << stepsize << endl;
  sal << "Modules: " << modules << endl;
  sal << endl << endl;
}


int Cap_evolI::make_modular_nw_unos_mues(GraphI &red, int *s, int nodes)
{
  if(!yapart_mues){
    cout << "[Error]: Partition was not initialized when  Cap_evolI::make_modular_nw_mues was called.\n";
    exit(1);
  }
    
  mu_mgood = new int [nodes];
  mu_mbad = new int [nodes];
  basic.create_array(mu_oc, nodes, modules);
  basic.create_array(mu_ot, nodes, modules);
  set<int>::iterator ite;
  set<int> whi;
  
  int count1 = 0;
  for(int n = 0; n < nodes; n++){ 
    mu_mbad[n] = 0;
    mu_mgood[n] = 0;
    for (int i = 0; i < modules; i++) {
        mu_oc[n][i] = 0;
        mu_ot[n][i] = 0;
        whi = *modkey[i];
        if(whi.count(n)!=0){
            for (ite = whi.begin(); ite != whi.end(); ite++){
                if(*ite >= nodes) {
                    cout << "[Error]: nodes in partition does not correspond to nodes in network in Cap_evolI::make_modular_nw.\n";
                    exit(1);
                }
            if(s[n]*s[*ite] == 1)
                red.change_interaction(n, *ite, 1);
            else
                red.change_interaction(n, *ite, -1);
            count1++;
            mu_mgood[n]++;
            }
        }
        whi.clear();
    }
  }
  mues_yaint = true;
  return count1;
}

void Cap_evolI::complete_conectiviness(GraphI& vacia, int *s, int nodes, int& nu_unos)
{
    int lugar_i, lugar_j;
    int valor;
    int ed;
    bool law;
    int chan_mod1;
    int chan_mod2;
    int mod_bu = 0;
    int mod_ma = 0;
    int *ocp;
    ocp = new int[modules];
    int *otp;
    otp = new int[modules];
    ed = vacia.number_of_edges();
    GraphI new_one(nodes, est, true);
    FitnessI lafi(est);
    new_one.copy(vacia);
    set<int> whi;
    
    if (ed < conectividad){
        for (int i = 0; i < conectividad - ed; i++){
             do{
                new_one.recorre_0((nodes*nodes - nu_unos), lugar_i, lugar_j);
                if (est.toss()) valor = 1;
                else valor = -1;
    
                new_one.change_interaction(lugar_j, lugar_i, valor);
                
                for (int k = 0; k < modules; k++) {
                    whi = *modkey[k];
                    if(whi.count(lugar_i)!=0)
                        chan_mod1 = k;
                    if(whi.count(lugar_j)!=0)
                        chan_mod2 = k;
                }
                
                mod_bu = mu_mgood[lugar_i];
                mod_ma = mu_mbad[lugar_i];
                if(chan_mod1 == chan_mod2){
                    if(s[lugar_i]*s[lugar_j]*valor == 1)
                        mod_bu++;
                    else 
                        mod_ma++;
                }
                
                for (int j = 0; j < modules; j++) {
                    ocp[j] = mu_oc[lugar_i][j];
                    otp[j] = mu_ot[lugar_i][j];
                }
                    
                if(chan_mod1 != chan_mod2){
                    if(s[lugar_i]*s[lugar_j]*valor == 1)
                        ocp[chan_mod2]++;
                    else 
                        otp[chan_mod2]++;
                }
                    
                law = lafi.cap_to_s_1gen(lugar_i, mod_bu, mod_ma, ocp, otp, modules);

                 if (law){
                     vacia.clear();
                     vacia.copy(new_one);
                     nu_unos++;
                     mu_mgood[lugar_i] = mod_bu;
                     mu_mbad[lugar_i] = mod_ma;
                     mu_oc[lugar_i][chan_mod2] = ocp[chan_mod2];
                     mu_ot[lugar_i][chan_mod2] = otp[chan_mod2];
                 }

                 else 
                    new_one.change_interaction(lugar_j, lugar_i, 0);
        
            } while (!law);
        }
    }
    
    else if (ed > conectividad){
        for (int i = 0; i < ed - conectividad; i++){
            do{
                new_one.recorre_1(nu_unos, lugar_i, lugar_j);
                valor = new_one.weight(lugar_j, lugar_i);

                new_one.change_interaction(lugar_j, lugar_i, 0);
 
                for (int k = 0; k < modules; k++) {
                    whi = *modkey[k];
                    if(whi.count(lugar_i)!=0)
                        chan_mod1 = k;
                    if(whi.count(lugar_j)!=0)
                        chan_mod2 = k;
                }
                
                mod_bu = mu_mgood[lugar_i];
                mod_ma = mu_mbad[lugar_i];
                if(chan_mod1 == chan_mod2){
                    if(s[lugar_i]*s[lugar_j]*valor == 1)
                        mod_bu --;
                    else 
                        mod_ma --;
                }
                
                for (int j = 0; j < modules; j++) {
                    ocp[j] = mu_oc[lugar_i][j];
                    otp[j] = mu_ot[lugar_i][j];
                }
                    
                if(chan_mod1 != chan_mod2){
                    if(s[lugar_i]*s[lugar_j]*valor == 1)
                        ocp[chan_mod2]--;
                    else 
                        otp[chan_mod2]--;
                }
                    
                law = lafi.cap_to_s_1gen(lugar_i, mod_bu, mod_ma, ocp, otp, modules);

                if (law){
                     vacia.clear();
                     vacia.copy(new_one);
                     nu_unos--;
                     mu_mgood[lugar_i] = mod_bu;
                     mu_mbad[lugar_i] = mod_ma;
                     mu_oc[lugar_i][chan_mod2] = ocp[chan_mod2];
                     mu_ot[lugar_i][chan_mod2] = otp[chan_mod2];
                 }

                 else 
                    new_one.change_interaction(lugar_j, lugar_i, valor);
                 
            } while (!law);
        }
    }
    
    whi.clear();
    new_one.clear();
    delete[] ocp;
    delete[] otp;
}

void Cap_evolI::muestreo_rn_v2_onenw(GraphI &vacia, int *s, int nodes, int &nu_unos, GraphI &regreso, int* old_mgood, int* old_mbad, int** old_oc, int** old_ot)
{  
//     cout << nu_unos << endl;
    int path = stepsize;
    bool law;
    int nu_unos_p;
    GraphI new_one(est);
    FitnessI lafi(est);
    
    int lugar_i, lugar_j, oldval, newval, chan_mod1, chan_mod2;
    int mod_bu = 0;
    int mod_ma = 0;
    int *ocp;
    ocp = new int[modules];
    int *otp;
    otp = new int[modules];
    
    nu_unos_p = nu_unos;
    new_one.copy(vacia);
    set<int> whi;
    
//         cout << conectividad << endl;

     for (int i = 1; i <= path; i++)
     {
         vacia.mutate_1g_v2(conectividad, nu_unos_p, lugar_i, lugar_j, oldval, newval);

         for (int k = 0; k < modules; k++) {
              whi = *modkey[k];
              if(whi.count(lugar_i)!=0)
                chan_mod1 = k;
              if(whi.count(lugar_j)!=0)
                chan_mod2 = k;
         }
         mod_bu = old_mgood[lugar_i];
         mod_ma = old_mbad[lugar_i];
         
         for (int j = 0; j < modules; j++) {
             ocp[j] = old_oc[lugar_i][j];
             otp[j] = old_ot[lugar_i][j];
         }
         
         if(oldval == 0){
            if(chan_mod1 == chan_mod2){
                if(s[lugar_i]*s[lugar_j]*newval == 1)
                    mod_bu++;
                else 
                    mod_ma++;
            }
            else{
                    if(s[lugar_i]*s[lugar_j]*newval == 1)
                        ocp[chan_mod2]++;
                    else 
                        otp[chan_mod2]++;
            }
         }
         
         if(oldval != 0){
             if(chan_mod1 == chan_mod2){
                if(s[lugar_i]*s[lugar_j]*oldval == 1)
                    mod_bu--;
                else 
                    mod_ma--;
            }
            else{
                    if(s[lugar_i]*s[lugar_j]*oldval == 1)
                        ocp[chan_mod2]--;
                    else 
                        otp[chan_mod2]--;
            }
         }
         
         law = lafi.cap_to_s_1gen(lugar_i, mod_bu, mod_ma, ocp, otp, modules);
        
         if (law){
//              cout << "si" << endl;
             new_one.clear();
             new_one.copy(vacia);
             nu_unos = nu_unos_p;
             old_mgood[lugar_i] = mod_bu;
             old_mbad[lugar_i] = mod_ma;
             old_oc[lugar_i][chan_mod2] = ocp[chan_mod2];
             old_ot[lugar_i][chan_mod2] = otp[chan_mod2];
          }
          else {
//               cout << "si" << endl;
              vacia.clear();
              vacia.copy(new_one);
              nu_unos_p = nu_unos;
          }
     }
    
    regreso.copy(new_one);
    new_one.clear();
    whi.clear();
    delete[] ocp;
    delete[] otp;
}

void Cap_evolI::muestreo_rn_v2_onenw_gaps(GraphI &vacia, int *s, int nodes, int &nu_unos, GraphI &regreso, int* old_mgood, int* old_mbad, int** old_oc, int** old_ot, int ***intmods, int nugaps)
{  
//     cout << nu_unos << endl;
    int path = stepsize;
    bool law;
    int nu_unos_p;
    GraphI new_one(est);
    FitnessI lafi (est);
    
    int lugar_i, lugar_j, oldval, newval, chan_mod1, chan_mod2;
    int mod_bu = 0;
    int mod_ma = 0;
    int *ocp;
    ocp = new int[modules];
    int *otp;
    otp = new int[modules];
    
    nu_unos_p = nu_unos;
    new_one.copy(vacia);
    set<int> whi;
    
//         cout << conectividad << endl;

     for (int i = 1; i <= path; i++)
     {
         vacia.mutate_1g_v2(conectividad, nu_unos_p, lugar_i, lugar_j, oldval, newval);

         for (int k = 0; k < modules; k++) {
              whi = *modkey[k];
              if(whi.count(lugar_i)!=0)
                chan_mod1 = k;
              if(whi.count(lugar_j)!=0)
                chan_mod2 = k;
         }
         mod_bu = old_mgood[lugar_i];
         mod_ma = old_mbad[lugar_i];
         
         for (int j = 0; j < modules; j++) {
             ocp[j] = old_oc[lugar_i][j];
             otp[j] = old_ot[lugar_i][j];
         }
         
         if(oldval == 0){
            if(chan_mod1 == chan_mod2){
                if(s[lugar_i]*s[lugar_j]*newval == 1)
                    mod_bu++;
                else 
                    mod_ma++;
            }
            else{
                    if(s[lugar_i]*s[lugar_j]*newval == 1)
                        ocp[chan_mod2]++;
                    else 
                        otp[chan_mod2]++;
            }
         }
         
         if(oldval != 0){
             if(chan_mod1 == chan_mod2){
                if(s[lugar_i]*s[lugar_j]*oldval == 1)
                    mod_bu--;
                else 
                    mod_ma--;
            }
            else{
                    if(s[lugar_i]*s[lugar_j]*oldval == 1)
                        ocp[chan_mod2]--;
                    else 
                        otp[chan_mod2]--;
            }
         }
         
         law = lafi.cap_to_s_1gen_gaps(lugar_i, mod_bu, mod_ma, ocp, otp, modules, intmods[lugar_i], nugaps);
        
         if (law){
//              cout << "si" << endl;
             new_one.clear();
             new_one.copy(vacia);
             nu_unos = nu_unos_p;
             old_mgood[lugar_i] = mod_bu;
             old_mbad[lugar_i] = mod_ma;
             old_oc[lugar_i][chan_mod2] = ocp[chan_mod2];
             old_ot[lugar_i][chan_mod2] = otp[chan_mod2];
          }
          else {
//               cout << "si" << endl;
              vacia.clear();
              vacia.copy(new_one);
              nu_unos_p = nu_unos;
          }
     }
    
    regreso.copy(new_one);
    new_one.clear();
    whi.clear();
    delete[] ocp;
    delete[] otp;
}

void Cap_evolI::start_gaps_ics(int **gaps, int nugaps, int modules, int nodes, int** modgapsp)
{
    
  int **modgaps;
  basic.create_array(modgaps, nugaps, modules);
  
  for(int i = 0; i < nugaps; i++){
      for(int j = 0; j < modules; j++){
        if(i == 0)
            modgaps[i][j] = 0;
        else
            modgaps[i][j] = abs(modgapsp[i][j] - modgapsp[0][j]); 
      }
  }
      
  int** s;
  basic.create_array(s, 2, nodes);
  for (int i = 0; i < nodes; i++) {
      if ((i%2)==0)
        s[0][i] = 1;
      else
        s[0][i] = -1;
      s[1][i] = s[0][i]*(-1);
    }
    
    
    int ***ms;
    basic.create_array(ms, 2, modules, nodes/modules);
    
    for(int i = 0; i < modules; i++)
        for(int j = 0; j < nodes/modules; j++){
            ms[0][i][j] = s[0][i*(nodes/modules) + j];
            ms[1][i][j] = s[1][i*(nodes/modules) + j];
        }
    
    int count;
    for(int j = 0; j < nugaps; j++)
        for (int i = 0; i < modules; i++) {
            count = 0;
            if (modgaps[j][i] == 1){
                for(int k = i*(nodes/modules); k < (i + 1)*(nodes/modules); k++){
                    gaps[j][k] = ms[1][i][count];
                    count++;
                }
            }
            else
                for(int k = i*(nodes/modules); k < (i + 1)*(nodes/modules); k++){
                    gaps[j][k] = ms[0][i][count];
                    count++;
                }
    }
}

void Cap_evolI::start_all_ics(int **alls, int nuic, int modules, int nodes)
{
  int** s;
  basic.create_array(s, 2, nodes);
  for (int i = 0; i < nodes; i++) {
      if ((i%2)==0)
        s[0][i] = 1;
      else
        s[0][i] = -1;
      s[1][i] = s[0][i]*(-1);
    }
    
    
    int ***ms;
    basic.create_array(ms, 2, modules, nodes/modules);
    
    for(int i = 0; i < modules; i++)
        for(int j = 0; j < nodes/modules; j++){
            ms[0][i][j] = s[0][i*(nodes/modules) + j];
            ms[1][i][j] = s[1][i*(nodes/modules) + j];
        }
    
    string *sbin;
    sbin = new string[nuic];
    for(int i = 0; i < nuic; i++){
        sbin[i] = basic.bintostring(i, modules);
    }
    
    
    int count;
    for(int j = 0; j < nuic; j++)
        for (int i = 0; i < modules; i++) {
            count = 0;
            if (sbin[j][i] == '1'){
                for(int k = i*(nodes/modules); k < (i + 1)*(nodes/modules); k++){
                    alls[j][k] = ms[1][i][count];
                    count++;
                }
            }
            else
                for(int k = i*(nodes/modules); k < (i + 1)*(nodes/modules); k++){
                    alls[j][k] = ms[0][i][count];
                    count++;
                }
    }
}
