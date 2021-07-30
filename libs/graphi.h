#ifndef GRAPHI_H
#define GRAPHI_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"


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
using std::setw;


class GraphI
{
public:
GraphI();
GraphI(Alea& jacta);
void start_rng(Alea& jacta);
GraphI(int n, Alea& jacta, bool dirnw);

//make network
void make_nw(int n, bool dirnw);
void copy(GraphI &templ);
void make_copy(GraphI& vacia);

//clear
void clear();
void clear_adj();
void clear_deg();
void clear_degdist();
void clear_attractor();
void clear_attractors();

//network's data
int number_of_edges();
void get_all_degrees();
int calc_indegree(int no);
int calc_outdegree(int no);
int calc_degree(int no);
int number_of_nodes();
int weight(int source, int target);
bool is_directed();

//analysis with networks
bool equal_nw(GraphI &templ);
double distance_between_nw(GraphI &templ);

//mutations
void recorre_0(int total, int& lugar_i, int& lugar_j);
void recorre_1(int total, int& lugar_i, int& lugar_j);
void change_interaction(int source, int target, int value);
void force_interaction(int source, int target, int value);
void mutate_1g_v2(int conectividad, int& nu_unos, int& lugar_i, int& lugar_j, int& oldval, int& newval);
void mutate(double con);
void mutate_gene(int gene, double con); 
void mutate_forced(double con);
void consider_mutation(double muratepg, double con);
void mate(GraphI &mother, GraphI &father);

//modularity
void build_moma_d(double **momadi);
void get_adjacency_matrix();
void copy_adjacency_matrix(int **matvac);

double eval_mod(double **mamod, set<set<int> > &equipos);
void partition_to_vector(const set<set<int> > &equipos, int *vop);
double eval_mod(double **mamod, int *vop);

//attractors and ability to sustain attractors
void interactions_to_attractor(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar);
void interactions_to_attractor_1gen(int *s, int mgood, int mbad, int* oc, int* ot, int nodes, int modules, set<set<int> >& predpar, int cu_gen);

void set_as_state(int *vec);
bool find_an_attractor_trunc(int steps);
void synchrony();

int attractor_size();
int attractor_element(int row, int node);

double distance_pub(int *goal);
double distance_pub(int **goal, int nugs);
double distance_aux_pub(int **goal, int nugs, int **atr, int numatr, int tam);

//read-write networks
void get_dir_nw_from_file(int nn, string arch);
void get_dir_nw_from_file(int nn, istream& en);
void get_dir_nw_from_file_bignw(int nn, istream& en, ifstream& index, int posicion, int elementos);
unsigned long get_address_from_index(ifstream& inventario, int posicion, int elementos);
void printnw(ostream& sal);
void export_nw(string arch);
void export_nw(ostream& fs);
void export_nw_fixed_width_bignw(ostream& fs, ofstream& key, int element);
void print_dot(string archi);
void print_dot_circ(string archi, double radius);
//print attractor
void print_attractor(ostream& sal);

//phenoypic variability
bool robustness_1mutation_gaps(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int nugaps, double &rob);
bool robustness_1mutation_cap(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, double &rob);

bool rob_by_dist_1mutation_gaps(int**ats, int nugaps, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double &rob, double &dist, double &robprom, double &robgan, double &robper, double &dist_gan, double &dist_per, double &mrob_gan, double &mrob_per, int &count_gan, int &count_per, int &trunc);
bool rob_by_dist_1mutation_cap(int**ats, int nic, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double &rob, double &dist, double &robprom,double &robgan, double &robper, double &dist_gan, double &dist_per, double &mrob_gan, double &mrob_per, int &count_gan, int &count_per, int &trunc);


void access_phen_1mutation_gaps(int**ats, int nugaps, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double prop_max_ac, int *tot_pht, double *simpht,  double *modspht, int *pht_g1, int *pht_d1, int *pht_10, int *pht_m1, int *pht_multchan, int *tot_acph, double *simacph, double *modsacph, int *acph_g1, int *acph_d1, int *acph_10, int *acph_m1, int *acph_multchan, int &maxsize, double &phen_size, double &dist_acph, int &trunc);
void access_phen_1mutation_cap(int**ats, int nuic, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double prop_max_ac, int *tot_pht, double *simpht,  double *modspht, int *pht_g1, int *pht_d1, int *pht_m1, int *pht_multchan, int *tot_acph, double *simacph, double *modsacph, int *acph_g1, int *acph_d1, int *acph_m1, int *acph_multchan, int &maxsize, double &phen_size, double &dist_acph, int &trunc);

void access_phen_translape_gaps(GraphI &temp, int**ats, int nugaps, int *s, int** mgood, int** mbad, int*** oc, int*** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double prop_max_ac, int **tot_phen, int **tot_acphen, int *tot_acphens, double *sim_acphens,  double *mods_acphens, int *acphens_d1, int *acphens_m1, int *acphens_multchan, int *tot_newpht, double *simnewpht,  double *modsnewpht, int *newpht_g1, int *newpht_d1, int *newpht_m1, int *newpht_multchan, int *tot_newacph, double *simnewacph, double *modsnewacph, int *newacph_g1, int *newacph_d1, int *newacph_m1, int *newacph_multchan, int &trunc);
void access_phen_translape_cap(GraphI &temp, int**ats, int nuic, int *s, int** mgood, int** mbad, int*** oc, int*** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double prop_max_ac, int **tot_phen, int **tot_acphen, int *tot_acphens, double *sim_acphens,  double *mods_acphens, int *acphens_d1, int *acphens_m1, int *acphens_multchan, int *tot_newpht, double *simnewpht,  double *modsnewpht, int *newpht_g1, int *newpht_d1, int *newpht_m1, int *newpht_multchan, int *tot_newacph, double *simnewacph, double *modsnewacph, int *newacph_g1, int *newacph_d1, int *newacph_m1, int *newacph_multchan, int &trunc);

int how_many_perturbed_modules(int *attorig, int modules, int& total_chan, double &sim);

private:
void set_default_exclusive_vars();
void prepare_for_degrees();
void prepare_for_degdist();


Alea est;
Basics basic;
int size;
int **nw;


bool directed;

int **matadya;
bool yamadya;

int numofe;
bool yae;
int *outdegree;
int *indegree;
int *degree;
bool yadeg;
int *degdist;
int *odegdist;
int *idegdist;
bool yadegdist;
    
int *edo;
int *ima;
int **attractor;
int atsize;
int palen;
bool yatra;
int ***attractors;
int *atsizes;
int numatrs;
bool yatras;
};
#endif
