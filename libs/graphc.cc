#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphc.h"
#include "graphi.h"
#include <string>

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


GraphC::GraphC()
{
}

GraphC::GraphC(Alea& jacta)
{
  start_rng(jacta);
}

void GraphC::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

GraphC::GraphC(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphC::make_nw(int n, bool dirnw)
{
  int i;
  size = n;
  nw = new double*[size];
  for (i = 0; i < size; i++)
		nw[i] = new double[size];
  basic.fillmat0(nw, size, size);
  directed = dirnw;
  set_default_exclusive_vars();
}

void GraphC::copy(GraphC &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

void GraphC::copy(GraphI &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

void GraphC::make_copy(GraphC& vacia)
{
  vacia.make_nw(size, is_directed());
  int i,j;
	if (is_directed()) {
		for (i=0; i<size; i++)
			for (j=0; j<size; j++)
				vacia.force_interaction(j, i, weight(j, i));
	}
	else {
		cout << "[Error]: GraphC::make_copy does not work with undirected graphs.\n";
		exit(1);
	}
}

void GraphC::clear()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] nw[i];
  delete [] nw;

  if (yadeg)
    clear_deg();
  if (yadegdist)
    clear_degdist();
  if (yamadya)
    clear_adj();
}

void GraphC::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphC::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
  if (yadegdist)
    clear_degdist();
}

void GraphC::clear_degdist()
{
	delete [] degdist;
	delete [] odegdist;
	delete [] idegdist;
	yadegdist = false;
}


int GraphC::number_of_edges()
{
  if (!yae) {
		numofe = 0;
		int i, j;
		if (directed) {
			for (i = 0; i < size; i++)
				for (j = 0; j < size; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		else {
			for (i = 0; i < size; i++)
				for (j = 0; j <= i; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		yae = true;
	}
  return numofe;
}

void GraphC::get_all_degrees()
{
  if (yadeg) {
		cout << "[Error]: Attempting to create degree arrays again. GraphC::get_all_degrees.\n";
		exit(1);
	}
  int i;
  prepare_for_degrees();
  for (i = 0; i < size; i++) {
		indegree[i] = calc_indegree(i);
		outdegree[i] = calc_outdegree(i);
		degree[i] = calc_degree(i);
	}
  yadeg = true;
}


int GraphC::calc_indegree(int no)
{
  int id = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[no][i] != 0)
      id++;
  return id;
}

int GraphC::calc_outdegree(int no)
{
  int od = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[i][no] != 0)
      od++;
  return od;
}

int GraphC::calc_degree(int no)
{
	int d;
	if (directed)
		d = indegree[no] + outdegree[no];
	else {
		if (indegree[no] != outdegree[no]) {
			cout << "[Error]: Different indegree and outdegree in an undirected network. GraphC::calc_degree.\n";
			exit(1);
		}
		d = indegree[no];
		if (weight(no, no) != 0)
			d++;
	}
	return d;
}

int GraphC::number_of_nodes()
{
  return size;
}

double GraphC::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphC::is_directed()
{
  return directed;
}

//analysis between networks
bool GraphC::equal_nw(GraphC &templ)
{
  int i, j;
  bool res = true;
  if ((templ.number_of_nodes() != size) || (templ.is_directed() != is_directed()))
    res = false;
  else {
    for (i = 0; i < size; i++) {
			for (j=0; j<size; j++)
				if (templ.weight(j, i) != weight(j, i)) {
					res = false;
					break;
				}
			if (!res)
				break;
		}
  }
  return res;
}

double GraphC::distance_between_nw(GraphC &templ)
{
  int i, j;
  double equals = 0;
  if ((templ.number_of_nodes() != size) || (templ.is_directed() != is_directed()))
    {
    cout << "[Error]: GraphI::compare_nw only compares netwoks of same size and only nw which are both directed or undirected\n";
    exit(1); 
    }
  else
    {
    for (i=0; i<size; i++)
	  for (j=0; j<size; j++)
	    equals+= abs(templ.weight(i, j) - weight(i, j));
    }
  return equals;
}

double GraphC::distance_between_nw(GraphI &templ)
{
  int i, j;
  double equals = 0;
  if ((templ.number_of_nodes() != size) || (templ.is_directed() != is_directed()))
    {
    cout << "[Error]: GraphI::compare_nw only compares netwoks of same size and only nw which are both directed or undirected\n";
    exit(1); 
    }
  else
    {
    for (i=0; i<size; i++)
	  for (j=0; j<size; j++)
	    equals+= abs(templ.weight(i, j) - weight(i, j));
    }
  return equals;
}

//modification
void GraphC::change_interaction(int source, int target, double value)
{
	if (directed) {
		if (value == nw[target][source]) {
			cout << "[Error]: GraphC::change_interaction changes to stay the same.\n";
			exit(1);
		}
		if (yae) {
      if ((value == 0) && (nw[target][source] != 0))
				numofe--;
      else
				if ((value != 0) && (nw[target][source] == 0))
					numofe++;
		}
		if (yadeg) {
      if ((value == 0) && (nw[target][source] != 0)) {
				outdegree[source]--;
				indegree[target]--;
				degree[source]--;
				degree[target]--;
			}
      else {
				if ((value != 0) && (nw[target][source] == 0)) {
					outdegree[source]++;
					indegree[target]++;
					degree[source]++;
					degree[target]++;
				}
      }
		}
		nw[target][source] = value;
	}
	else {
		cout << "[Error]: GraphC::change_interaction does not work for undirected graphs.\n";
		exit(1);
  }
  if (yadegdist)
    clear_degdist();
  if (yamadya)
    clear_adj();
}

void GraphC::force_interaction(int source, int target, double value)
{
	if (directed) {
		if (yae) {
      if ((value == 0) && (nw[target][source] != 0))
				numofe--;
      else
				if ((value != 0) && (nw[target][source] == 0))
					numofe++;
		}
		if (yadeg) {
      if ((value == 0) && (nw[target][source] != 0)) {
				outdegree[source]--;
				indegree[target]--;
				degree[source]--;
				degree[target]--;
			}
      else {
				if ((value != 0) && (nw[target][source] == 0)) {
					outdegree[source]++;
					indegree[target]++;
					degree[source]++;
					degree[target]++;
				}
      }
		}
		nw[target][source] = value;
	}
	else {
		cout << "[Error]: force_interaction does not work for undirected graphs. GraphC::force_interaction.\n";
		exit(1);
  }
  if (yadegdist)
    clear_degdist();
  if (yamadya)
    clear_adj();
}

//private
void GraphC::set_default_exclusive_vars()
{
  yae = false;
  yadeg = false;
  yadegdist=false;
  yamadya = false;
}

void GraphC::prepare_for_degrees()
{
	degree = new int[size];
	outdegree = new int[size];
	indegree = new int[size];
}

void GraphC::prepare_for_degdist()
{
	degdist = new int[size+1];
	basic.fillv0(degdist, size+1);
	
	odegdist = new int[size+1];
	basic.fillv0(odegdist, size+1);
	
	idegdist = new int[size+1];
	basic.fillv0(idegdist, size+1);
}
