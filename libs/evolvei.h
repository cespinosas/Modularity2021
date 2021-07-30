#ifndef EVOLVEI_H
#define EVOLVEI_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi.h"
#include "graphc.h"
#include "fitnessi.h"

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

class EvolveI
{
public:
	EvolveI();
  EvolveI(Alea& jacta);
  void start_rng(Alea& jacta);
  void set_params(double muratep, double wolp, int popsizep, int nugenp);
  void print_params(ostream& sal);
  void start_pop(GraphI &founder);
  
  void one_generation();
  void one_generation_mod(double modmax, double modmin, set<set<int> >& predparp);
  int one_generation_mod_dead(double modmax, double modmin, set<set<int> >& predparp);
  void one_generation_sex();
  void assign_w(int cual, double cuanto);
  void calc_meanw();
  double return_meanw();
  double return_maxw();
  double return_sdw();
  void optima(set<int> &opt);
  void optima_strict(set<int> &opt);
  int num_optima();
  int num_optima_strict(); 
  void clear();
  double get_murate();
  double get_wol();
  int get_popsize();
  int get_nugen();
  double get_w(int cual);
  
  void mean_nw(GraphC &nw);
  double more_repeated_nw(GraphI & nw);
  
  double dist_bw_nws();
  
  GraphI *population;
  
private:

  Alea est;
  Basics basic;
  double murate;
  double wol;
  int popsize;
  int nugen;
  
  bool yamean;
  double meanw;
  double maxw;
  double sdw;
  double *w;
  GraphI *nepop;
  bool yapop;
  bool yaparams;
  
};

#endif
