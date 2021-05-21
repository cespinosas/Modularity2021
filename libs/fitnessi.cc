#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "fitnessi.h"
#include <string>
#include "graphi.h"

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


FitnessI::FitnessI()
{
}

FitnessI::FitnessI(Alea& jacta)
{
  start_rng(jacta);
}

void FitnessI::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

bool FitnessI::cap_to_s(int nodes, int* mod_bu, int* mod_ma, int **oc, int **ot, int  modules)
{
    bool flag = true;
    for(int n = 0; n < nodes; n++){
       flag = cap_to_s_1gen(n, mod_bu[n], mod_ma[n], oc[n], ot[n], modules);
       if (!flag)
           break;
    }
    return flag;
}

bool FitnessI::cap_to_s_1gen(int gen, int mod_bu, int mod_ma, int *oc, int *ot, int  modules)
{
    int q = 0;
    for (int j = 0; j < modules; j++)
        q += abs(oc[j]-ot[j]);
        
    if(mod_bu >= (mod_ma + q))
        return true;
    else
        return false;
    
}

bool FitnessI::cap_to_s_gaps(int nodes, int* mod_bu, int* mod_ma, int **oc, int **ot, int  modules, int ***intmods, int nugaps)
{
    bool flag = true;
    for(int n = 0; n < nodes; n++){
       flag = cap_to_s_1gen_gaps(n, mod_bu[n], mod_ma[n], oc[n], ot[n], modules, intmods[n], nugaps);
       if (!flag)
           break;
    }
    return flag;
}


bool FitnessI::cap_to_s_1gen_gaps(int gen, int mod_bu, int mod_ma, int *oc, int *ot, int  modules, int **intmods, int nugaps)
{
    int q;
    int Nt = 0;
    for (int i = 0; i < nugaps; i++){
        q = 0;
        Nt = mod_bu - mod_ma;
        for (int j = 0; j < modules; j++)
            q = q + intmods[i][j]*(oc[j]-ot[j]);
        if(!(Nt + q >= 0))
            return false;
    }
        
    return true;
}

double FitnessI::distance(int *goal, int *atr, int tam)
{
  int dif = 0,i;
  for (i=0; i<tam; i++)
    if (goal[i] != atr[i])
      dif++;
  return double(dif)/double(tam);
}

double FitnessI::distance(int *goal, int **atr, int numatr, int tam)
{
  int i,j;
  double dif = 0;
  for (i=0; i<numatr; i++)
    for (j=0; j < tam; j++)
      if (goal[j] != atr[i][j])
        dif = dif + (1.0/double(numatr));
  return dif/double(tam);
}

double FitnessI::distance(int **goal, int nugs, int **atr, int numatr, int tam)
{
  int i,j,k;
  double dis;
  if (nugs != numatr) {
    int mulng, mulna;
    if (nugs < numatr) {
      if ((numatr%nugs)==0) {
        mulna = 1;
        mulng = numatr/nugs;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    else {
      if ((nugs%numatr)==0) {
        mulng = 1;
        mulna = nugs/numatr;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    int **nugoal, **nuatr;
    nugoal = new int*[nugs*mulng];
    for (i= 0; i<(nugs*mulng); i++)
      nugoal[i] = new int[tam];
    nuatr = new int*[numatr*mulna];
    for (i=0; i <(numatr*mulna); i++)
      nuatr[i] = new int[tam];
    for (i=0; i<mulng; i++)
      for (j=0; j<nugs; j++)
        for (k=0; k < tam; k++)
          nugoal[(nugs*i)+j][k] = goal[j][k];
    for (i=0; i < mulna; i++)
      for (j=0; j<numatr; j++)
        for (k=0; k<tam;k++)
          nuatr[(numatr*i)+j][k] = atr[j][k];
    dis = distance_aux(nugoal, (nugs*mulng), nuatr, (numatr*mulna), tam);
    for (i= 0; i<(nugs*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(numatr*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
  }
  else
    dis = distance_aux(goal, nugs, atr, numatr, tam);
  return dis;
}

double FitnessI::distance_aux(int **goal, int nugs, int **atr, int numatr, int tam)
{
  if (nugs != numatr) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists,d;
  dists = new double[nugs];
  basic.fillv0(dists, nugs);
  for (i=0; i<nugs; i++) {
    dists[i] =0;
    for (j=0; j<nugs; j++) {
      for (k=0;k<tam;k++)
        if (goal[j][k] != atr[(j+i)%nugs][k])
          dists[i] = dists[i] + (1.0/double(nugs));
    }
  }
  d = basic.find_min(dists, nugs);
  d = d/double(tam);
  delete [] dists;
  return d;
}
