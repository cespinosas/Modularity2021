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
#include "graphc.h"
#include "fitnessi.h"
#include "evolvei.h"

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

EvolveI::EvolveI()
{
  yaparams = false;
  yapop = false;
}

EvolveI::EvolveI(Alea& jacta)
{
  start_rng(jacta);
  yaparams = false;
  yapop = false;
}

void EvolveI::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

void EvolveI::set_params(double muratep, double wolp, int popsizep, int nugenp)
{
  murate = muratep;
  wol = wolp;
  popsize = popsizep;
  nugen = nugenp;
  yaparams = true;
}

void EvolveI::print_params(ostream& sal)
{
  sal << "Integer weights\n";
  sal << "Number of genes: " << nugen << endl;
  sal << "Population size: " << popsize << endl;
  sal << "Mutation rate: " << murate << endl;
  sal << "Expected connectivity: " << wol << endl;
  sal << endl << endl;
}

void EvolveI::start_pop(GraphI &founder)
{
  if (!yaparams) {
    cout << "[Error]: Parameters had not been set when EvolveI::start_pop was called.\n";
    exit(1);
  }
  if (yapop) {
    cout << "[Error]: Population already constructed when EvolveI::start_pop was called.\n";
    exit(1);
  }
  population = new GraphI[popsize];
  nepop = new GraphI[popsize];
  w = new double[popsize];
  basic.fillv0(w, popsize);
  meanw =0;
  maxw = 0;
  int i;
  for (i=0; i<popsize; i++) {
    population[i].start_rng(est);
    nepop[i].start_rng(est);
    population[i].copy(founder);
    population[i].consider_mutation(murate, wol);
  }
  yapop = true;
  yamean = false;
}

void EvolveI::one_generation()
{
  int i,j;
  double frog, toad;
  for (i=0; i<popsize; i++) {
    if (!yamean){
        cout << "[Error]: meanw had not been calculated when EvolveI::one generation was called .\n";
          exit(1);
    }
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    nepop[i].copy(population[j]);
    nepop[i].consider_mutation(murate, wol);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
  yamean = false;
}

void EvolveI::one_generation_mod(double modmax, double modmin, set<set<int> >& predparp)
{
  int j;
  double frog, toad;
  double mod;
  double **matmod;
  basic.create_array(matmod, nugen, nugen);
  GraphI temp(est);
  bool fmod;
  for (int i=0; i<popsize; i++) {
    if (!yamean){
        cout << "[Error]: meanw had not been calculated when EvolveI::one generation was called .\n";
          exit(1);
    }
    
    fmod = false;
    do{
        frog = est.randreal()*meanw*popsize;
        toad = 0;
        for (j=0; j<popsize; j++) {
            toad = toad + w[j];
            if (toad >= frog)
                break;
        }
        temp.copy(population[j]);
        temp.consider_mutation(murate, wol);
        temp.build_moma_d(matmod);
        mod = temp.eval_mod(matmod, predparp);
        
        if(mod < modmax && mod > modmin){
            nepop[i].copy(temp);
            fmod = true;
        }
        temp.clear();
        
    }while(!fmod);
  }
  
  for (int i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
  yamean = false;
  for(int i = 0; i < nugen; i++)
      delete matmod[i];
  delete[] matmod;
}

int EvolveI::one_generation_mod_dead(double modmax, double modmin, set<set<int> >& predparp)
{
  int dead = 0;
  int j;
  double frog, toad;
  double mod;
  double **matmod;
  basic.create_array(matmod, nugen, nugen);
  GraphI temp(est);
  bool fmod;
  for (int i=0; i<popsize; i++) {
    if (!yamean){
        cout << "[Error]: meanw had not been calculated when EvolveI::one generation was called .\n";
          exit(1);
    }
    
    fmod = false;
    do{
        frog = est.randreal()*meanw*popsize;
        toad = 0;
        for (j=0; j<popsize; j++) {
            toad = toad + w[j];
            if (toad >= frog)
                break;
        }
        temp.copy(population[j]);
        temp.consider_mutation(murate, wol);
        temp.build_moma_d(matmod);
        mod = temp.eval_mod(matmod, predparp);
        
        if(mod < modmax && mod > modmin){
            nepop[i].copy(temp);
            fmod = true;
        }
        else
            dead++;
        temp.clear();
        
    }while(!fmod);
  }
  
  for (int i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
  yamean = false;
  for(int i = 0; i < nugen; i++)
      delete matmod[i];
  delete[] matmod;
  return dead;
  
}

void EvolveI::one_generation_sex()
{
  int i,j,k;
  double frog, toad, tadp;
  for (i=0; i<popsize; i++) {
    if (!yamean){
        cout << "[Error]: meanw had not been calculated when EvolveI::one generation was called .\n";
          exit(1);
    }
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    frog = est.randreal()*meanw*popsize;
    tadp = 0;
    for (k=0; k< popsize; k++) {
      tadp = tadp + w[k];
      if (tadp >= frog)
        break;
    }
    nepop[i].mate(population[j], population[k]);
    nepop[i].consider_mutation(murate, wol);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
  yamean = false;
}

void EvolveI::assign_w(int cual, double cuanto)
{
  if ((cual >= popsize) || (cual < 0)) {
    cout << "[Error]: Individual does not exist. EvolveI::assign_w.\n";
    exit(1);
  }
  w[cual] = cuanto;
}

void EvolveI::calc_meanw()
{
  meanw = basic.get_mean(w, popsize);
  maxw = basic.find_max(w, popsize);
  sdw = basic.get_pop_stddev(w, popsize);
  yamean = true;
}

double EvolveI::return_meanw()
{
  if (!yamean){
    cout << "[Error]: meanw had not been calculated when EvolveI::return_meanw was called .\n";
    exit(1);
  }
  return meanw;
}

double EvolveI::return_maxw()
{
  if (!yamean){
    cout << "[Error]: maxw had not been calculated when EvolveI::return_maxw was called .\n";
    exit(1);
  }
  return maxw;
}

double EvolveI::return_sdw()
{
  if (!yamean){
    cout << "[Error]: maxw had not been calculated when EvolveI::return_maxw was called .\n";
    exit(1);
  }
  return sdw;
}

void EvolveI::optima(set<int> &opt)
{
  opt.clear();
  int i;
  for (i=0; i < popsize; i++)
    if (w[i]==maxw)
      opt.insert(i);
}

void EvolveI::optima_strict(set<int> &opt)
{
  opt.clear();
  int i;
  for (i=0; i < popsize; i++)
    if (w[i]==1)
      opt.insert(i);
}

int EvolveI::num_optima()
{
  int i, res=0;
  for (i=0; i < popsize; i++)
    if (w[i]==maxw)
      res++;
  return res;
}

int EvolveI::num_optima_strict()
{
  int i, res=0;
  for (i=0; i < popsize; i++)
    if (w[i]==1)
      res++;
  return res;
}

void EvolveI::clear()
{
  int i;
  if (yapop) {
    for (i=0; i<popsize; i++)
      population[i].clear();
    delete [] population;
    delete [] nepop;
    delete [] w;
  }
  yapop = false;
  yaparams = false;
  yamean = false;
}


double EvolveI::get_murate()
{
  return murate;
}

double EvolveI::get_wol()
{
  return wol;
}

int EvolveI::get_popsize()
{
  return popsize;
}

int EvolveI::get_nugen()
{
  return nugen;
}

double EvolveI::get_w(int cual) {
  return w[cual];
}


void EvolveI::mean_nw(GraphC &nw)
{
    double interac;
    
    for (int i = 0; i < nugen; i++){
        for(int j = 0; j < nugen; j++){
            interac = 0;
            for(int k = 0; k < popsize; k++){
                interac = interac + population[k].weight(j,i);
            }
            if(interac != 0)
                nw.change_interaction(j, i, (double)interac/popsize);
        }
    }
}

double EvolveI::more_repeated_nw(GraphI & nw)
{
  int *rep_times;
  rep_times = new int[popsize];
  bool *comp_nw;
  comp_nw = new bool[popsize];
  basic.fillv0(rep_times, popsize);
  basic.fillv1(comp_nw, popsize);
  
  for(int i = 0; i < popsize; i++){
      if(comp_nw[i]){
          for(int j = i+1; j < popsize; j++){
              if(comp_nw[j])
                  if(population[i].equal_nw(population[j])){
                      rep_times[i]++;
                      comp_nw[j] = false;
                  }
          }
      }
  }
  
  int cual;
  cual = basic.find_max_index(rep_times, popsize);
  nw.copy(population[cual]);
  
  return((double)basic.find_max(rep_times, popsize)/popsize);
}

double EvolveI::dist_bw_nws()
{
    double dist_nws = 0;
    for(int i = 0; i < popsize; i++){
        for(int j = (i+1); j < popsize; j++){
            dist_nws += population[i].distance_between_nw(population[j]);
        }
    }
    
    return dist_nws/(double((popsize*(popsize-1)))/2.0);
        
}
