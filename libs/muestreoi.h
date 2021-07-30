#ifndef MUESTREOI_H
#define MUESTREOI_H

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

class MuestreoI
{
public:

MuestreoI();
MuestreoI(Alea& jacta);
void start_rng(Alea& jacta);

//sample functions
void set_params_muestreo(int samplesizep, int nugenp, int stepsizep, int conectividadp, int modulesp);
void set_partition(set<set<int> > &predparp);
void print_params_mues(ostream& sal);

int make_modular_nw_unos_mues(GraphI &red, int *s, int nodes);
void complete_conectiviness(GraphI& vacia, int *s, int nodes, int& nu_unos);

void muestreo_rn_v2_onenw(GraphI &vacia, int *s, int nodes, int &nu_unos, GraphI &regreso, int* old_mgood, int* old_mbad, int** old_oc, int** old_ot);
void muestreo_rn_v2_onenw_gaps(GraphI &vacia, int *s, int nodes, int &nu_unos, GraphI &regreso, int* old_mgood, int* old_mbad, int** old_oc, int** old_ot, int ***intmods, int nugaps);


void start_gaps_ics(int **gaps, int nugaps, int modules, int nodes, int** modgapsp);
void start_all_ics(int **alls, int nuic, int modules, int nodes);

private:
Alea est;
Basics basic;
FitnessI law; 

int samplesize;
int stepsize;
int modules;
int conectividad;
int nugen;
bool yaparamsmu;

bool mues_yaint;
int *mu_mgood;
int *mu_mbad;
int **mu_oc;    
int **mu_ot;
bool yapart;
bool yapart_mues;
set<set<int> > predpar;
set<set<int> >::iterator *modkey;
};

#endif
