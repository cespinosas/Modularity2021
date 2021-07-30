#ifndef GRAPHC_H
#define GRAPHC_H

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
using std::setw;


class GraphC
{
public:
  GraphC();
  GraphC(Alea& jacta);
  void start_rng(Alea& jacta);
  GraphC(int n, Alea& jacta, bool dirnw);
  
  //make network
  void make_nw(int n, bool dirnw);
  void copy(GraphC &templ);
  void copy(GraphI &templ);
  void make_copy(GraphC& vacia);
  
  //clear
  void clear();
  void clear_adj();
  void clear_deg();
  void clear_degdist();
  
  //network's data
  int number_of_edges();
  void get_all_degrees();
  void get_degdist();
  int calc_indegree(int no); 
  int calc_outdegree(int no); 
  int calc_degree(int no);
  int number_of_nodes();
  double weight(int source, int target);
  bool is_directed();
  
  //analysis
  bool equal_nw(GraphC &templ);
  double distance_between_nw(GraphC &templ);
  double distance_between_nw(GraphI &templ);
  
  //modification
  void change_interaction(int source, int target, double value);
  void force_interaction(int source, int target, double value);
  
	
private:
    
  void set_default_exclusive_vars();
  void prepare_for_degrees();
  void prepare_for_degdist();
    
  Alea est;
  Basics basic;
  int size;
  double **nw;
  
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
 

};


#endif
