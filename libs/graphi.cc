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
#include "graphi.h"
#include "fitnessi.h"
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


GraphI::GraphI()
{
}

GraphI::GraphI(Alea& jacta)
{
  start_rng(jacta);
}

void GraphI::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

GraphI::GraphI(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphI::make_nw(int n, bool dirnw)
{
  int i;
  size = n;
  nw = new int*[size];
  for (i = 0; i < size; i++) 
		nw[i] = new int[size];
  basic.fillmat0(nw, size, size);
  directed = dirnw;
  set_default_exclusive_vars();
  edo = new int[size];
  ima = new int[size];
}

void GraphI::copy(GraphI &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

void GraphI::make_copy(GraphI& vacia)
{
  vacia.make_nw(size, is_directed());
  int i,j;
	if (is_directed()) {
		for (i=0; i<size; i++)
			for (j=0; j<size; j++)
				vacia.force_interaction(j, i, weight(j, i));
	}
	else {
		cout << "[Error]: GraphI::make_copy does not work with undirected graphs.\n";
		exit(1);
	}
}

void GraphI::clear()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] nw[i];
  delete [] nw;
  delete [] ima;
  delete [] edo;

  if (yadeg)
    clear_deg();
  if (yadegdist)
    clear_degdist();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yamadya)
    clear_adj();
}

void GraphI::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphI::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
  if (yadegdist)
      clear_degdist();
}

void GraphI::clear_degdist()
{
	delete [] degdist;
	delete [] odegdist;
	delete [] idegdist;
	yadegdist = false;
}

void GraphI::clear_attractor()
{
  int i;
  for (i = 0; i < atsize; i++)
    delete [] attractor[i];
  delete [] attractor;
  atsize = 0;
  palen = 0;
  yatra = false;
}

void GraphI::clear_attractors()
{
  int i,j;
  for (i=0; i < numatrs; i++) {
    for (j=0; j< atsizes[i]; j++)
      delete [] attractors[i][j];
    delete [] attractors[i];
  }
  delete [] attractors;
  delete [] atsizes;
  numatrs = 0;
  yatras = false;
}

int GraphI::number_of_edges()
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

void GraphI::get_all_degrees()
{
  if (yadeg) {
		cout << "[Error]: Attempting to create degree arrays again. GraphI::get_all_degrees.\n";
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

int GraphI::calc_indegree(int no)
{
  int id = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[no][i] != 0)
      id++;
  return id;
}

int GraphI::calc_outdegree(int no)
{
  int od = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[i][no] != 0)
      od++;
  return od;
}

int GraphI::calc_degree(int no)
{
  int d;
  if (directed)
    d = indegree[no] + outdegree[no];
  else {
    if (indegree[no] != outdegree[no]) {
			cout << "[Error]: Different indegree and outdegree in an undirected network. GraphI::calc_degree.\n";
			exit(1);
		}
    d = indegree[no];
    if (weight(no, no) != 0)
        d++;
  }
  return d;
}

int GraphI::number_of_nodes()
{
  return size;
}

int GraphI::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphI::is_directed()
{
  return directed;
}

bool GraphI::equal_nw(GraphI &templ)
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

double GraphI::distance_between_nw(GraphI &templ)
{
  int i, j;
  int equals = 0;
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

void GraphI::recorre_0(int total, int& lugar_i, int& lugar_j)
{
    bool flag = false;
    int i, j;
    int lugar;
    int count = 0;
    lugar = est.randint(1, total + 1);
    for (i = 0; i < size; i ++){
        for (j = 0; j < size; j++){
            if (weight(j, i) == 0)
                count++;
            if (count == lugar){
                flag = true;
                break;
            }
        }
            if (flag == true) break;
    }
    
    lugar_i = i;
    lugar_j = j;
}

void GraphI::recorre_1(int total, int& lugar_i, int& lugar_j)
{
    bool flag = false;
    int i, j;
    int lugar;
    int count = 0;
    lugar = est.randint(1, total + 1);
    for (i = 0; i < size; i ++){
        for (j = 0; j < size; j++){
            if (weight(j, i) != 0)
                count++;
            if (count == lugar){
                flag = true;
                break;
            }
        }
            if (flag == true) break;
    }
    
    lugar_i = i;
    lugar_j = j;
}

void GraphI::change_interaction(int source, int target, int value)
{
	if (directed) {
		if (value == nw[target][source]) {
			cout << "[Error]: GraphI::change_interaction changes to stay the same.\n";
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
		cout << "[Error]: GraphI::change_interaction does not work for undirected graphs.\n";
		exit(1);
  }
  if (yadegdist)
    clear_degdist();
  if (yamadya)
    clear_adj();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
}

void GraphI::force_interaction(int source, int target, int value)
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
		cout << "[Error]: GraphI::force_interaction does not work for undirected graphs.\n";
		exit(1);
  }
  if (yadegdist)
    clear_degdist();
  if (yamadya)
    clear_adj();
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
}

void GraphI::mutate_1g_v2(int conectividad, int& nu_unos, int& lugar_i, int& lugar_j, int& oldval, int& newval) {
  if (!directed) {
    cout << "[Error]: GraphC::mutate_gene_v2 is not defined for undirected networks.\n";
    exit(1);
  }
  int ed = number_of_edges();;
  double valor;
  
  if ((ed < (conectividad - 2)) || (ed > (conectividad + 2))) {
    cout << "[Error]: GraphI::mutate_gene_v2 is operable only for networks with pre-fixed conectivity.\n";
    exit(1);
  }
  
  if ((ed > (conectividad - 2)) && (ed < (conectividad + 2))){
      if (est.toss()){
          recorre_0((size*size - nu_unos), lugar_i, lugar_j);
          oldval = 0;
          if (est.toss()) valor = 1;
          else valor = -1;
          
          change_interaction(lugar_j, lugar_i, valor);
          newval = valor;
          nu_unos++;
      }
      else {
          recorre_1(nu_unos, lugar_i, lugar_j);
          
          oldval = weight(lugar_j, lugar_i);
          change_interaction(lugar_j, lugar_i, 0);
          newval = 0;
          nu_unos--;
      }
  }
    
  else if (ed == (conectividad - 2)){
          recorre_0((size*size - nu_unos), lugar_i, lugar_j);
          oldval = 0;
          if (est.toss()) valor = 1;
          else valor = -1;
          
          change_interaction(lugar_j, lugar_i, valor);
          newval = valor;
          nu_unos++;
      }

  else if (ed == (conectividad + 2)){
          recorre_1(nu_unos, lugar_i, lugar_j);
          
          oldval = weight(lugar_j, lugar_i);
          change_interaction(lugar_j, lugar_i, 0);
          newval = 0;
          nu_unos--;
  }

}

void GraphI::mutate(double con) 
{
  int i = est.randint(0, size);
  mutate_gene(i, con);
}

void GraphI::mutate_gene(int gene, double con) 
{
  if (!directed) {
    cout << "[Error]: GraphI::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  bool sube = false;
  int j;
  j = est.randint(0, size);
  int oldval = weight(j,gene);
  if (est.randreal() < con)
    sube = true;
  if (sube) {
    if (oldval == 0) {
      if (est.toss())
        force_interaction(j,gene,1);
      else
        force_interaction(j,gene,-1);
    }
  }
  else {
    if (oldval != 0)
      force_interaction(j, gene, 0);
  }
}

void GraphI::mutate_forced(double con) 
{
  if (!directed) {
    cout << "[Error]: GraphI::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  bool sube = false;
  int i, j, oldval;
  bool mut = false;
  
  while(!mut){
    i = est.randint(0, size);
    j = est.randint(0, size);
    oldval = weight(j, i);
    
    if (est.randreal() < con)
        sube = true;
    if (sube) {
        if (oldval == 0) {
            if (est.toss())
                force_interaction(j, i,1);
            else
                force_interaction(j, i,-1);
            mut = true;
        }
    }
    else {
        if(oldval != 0){
            force_interaction(j, i, 0);
            mut = true;
        }
        
    }
  }
}

//Para la evolución
void GraphI::consider_mutation(double muratepg, double con)
{
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene(i, con);
}

void GraphI::mate(GraphI &mother, GraphI &father)
{
  if (mother.number_of_nodes() != father.number_of_nodes()) {
		cout << "[Error]: Parents with different genome size. GraphI::mate.\n";
		exit(1);
	}
	if (mother.is_directed() != father.is_directed()) {
		cout << "[Error]: One parent is a directed network while the other is an undirected network. GraphI::mate.\n";
		exit(1);
	}
  int i, j;
  bool paoma;
  int n = father.number_of_nodes();
  make_nw(n, father.is_directed());
	if (directed) {
		for (i = 0; i < size; i++) {
			paoma = est.toss();
			if (paoma) {
				for (j = 0; j < size; j++)
          force_interaction(j,i,father.weight(j,i));
			}
			else {
				for (j = 0; j < size; j++)
          force_interaction(j,i,mother.weight(j,i));
			}
		}
	}
	else {
		for (i=0; i<size; i++)
			for (j=0; j<=i; j++) {
				paoma = est.toss();
				if (paoma)
          force_interaction(j,i,father.weight(j,i));
				else
          force_interaction(j,i,mother.weight(j,i));
			}
	}
}

//Modularidad
void GraphI::build_moma_d(double **momadi) {
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphI::build_moma_d\n";
    exit(1);
  }
  if (!yamadya)
    get_adjacency_matrix();
  if (!yadeg)
    get_all_degrees();
  int i, j;
  for (i = 0; i < size; i++)
    for (j= 0; j < size; j++)
      momadi[i][j] = matadya[i][j] + matadya[j][i] - ((indegree[i]*outdegree[j]) + (indegree[j]*outdegree[i]))/(double(number_of_edges()));
}

void GraphI::get_adjacency_matrix() {
  if (yamadya) {
    cout << "[Error]: Adjacency matrix already created when GraphI::get_adjacency_matrix was called.\n";
    exit(1);
  }
  basic.create_array(matadya, size, size);
  int i,j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++) {
      if (nw[i][j] != 0)
        matadya[i][j] = 1;
      else
        matadya[i][j] = 0;
    }
  if (!directed) {
    for (i = 0; i < size; i++)
      if (matadya[i][i] != 0)
        matadya[i][i]++;
  }
  yamadya = true;
  return;
}

void GraphI::copy_adjacency_matrix(int **matvac)
{
  if (!yamadya) {
      get_adjacency_matrix();
  }
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
        matvac[i][j] = matadya[i][j];
}

double GraphI::eval_mod(double **mamod, set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Q = eval_mod(mamod, arr);
  delete [] arr;
  return Q;
}

void GraphI::partition_to_vector(const set<set<int> > &equipos, int *vop) {
  set<set<int> >::iterator it;
  set<int>::iterator ite;
  set<int> whi;
  basic.fillvm1(vop, size);
  int conc = 0;
  for (it = equipos.begin(); it != equipos.end(); it++) {
    whi = *it;
    for (ite = whi.begin(); ite != whi.end(); ite++)
      vop[*ite] = conc;
    conc++;
  }
}

double GraphI::eval_mod(double **mamod, int *vop) {
  int i, j;
  double Q = 0;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (vop[i] == vop[j])
        Q += mamod[i][j];
  Q /= (2.0*number_of_edges());
  return Q;
}

//get interactions to mantain an atractor s
void GraphI::interactions_to_attractor(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar)
{
  set<set<int> >::iterator *modkey;
  modkey = new set<set<int> >::iterator [predpar.size()];
  set<int>::iterator ite;
  set<int> whi;
  
  if(modules != (int)predpar.size()){
      cout << "[Error]: Number of modules given and modules in partition are not the same in GraphI::iteractions_to_attractor.\n";
      exit(1);
  }
  
  for(int i = 0; i < modules; i++){
      modkey[i] = predpar.begin();
      for(int j = 0; j < i; j++)
          modkey[i]++;
  }
  
  int aij;
  for(int n = 0; n < nodes; n++){
        mgood[n] = 0;
        mbad[n] = 0;
        for (int i = 0; i < modules; i++) {
            oc[n][i] = 0;
            ot[n][i] = 0;
            whi = *modkey[i];
            if(whi.count(n)!=0){
                for (ite = whi.begin(); ite != whi.end(); ite++){
                    if(*ite >= nodes){
                        cout << "[Error]: Nodes in partition does not belong to network in GraphI::iteractions_to_attractor.\n";
                        exit(1);
                    }
                    aij = weight(*ite, n);
                    if(aij != 0){
                        if(s[n]*s[*ite]*aij == 1)
                            mgood[n]++;
                        else
                            mbad[n]++;
                    }
                }
                    
            }
            else{
                for (ite = whi.begin(); ite != whi.end(); ite++){
                    if(*ite >= nodes){
                        cout << "[Error]: Nodes in partition does not belong to network in GraphI::iteractions_to_attractor.\n";
                        exit(1);
                    }
                    aij = weight(*ite, n);                                        /*cout << aij << " " << *ite << "\t";*/
                    if(aij != 0){
                        if(s[n]*s[*ite]*aij == 1)
                            oc[n][i]++;
                        else
                            ot[n][i]++;
                    }
                }
            }
        }

    }
    whi.clear();
    delete[] modkey;
}
 
void GraphI::interactions_to_attractor_1gen(int *s, int mgood, int mbad, int* oc, int* ot, int nodes, int modules, set<set<int> >& predpar, int cu_gen)
{
  set<set<int> >::iterator *modkey;
  modkey = new set<set<int> >::iterator [predpar.size()];
  set<int>::iterator ite;
  set<int> whi;
  
  if(modules != (int)predpar.size()){
      cout << "[Error]: Number of modules given and modules in partition are not the same in GraphI::iteractions_to_attractor.\n";
      exit(1);
  }
  
  for(int i = 0; i < modules; i++){
      modkey[i] = predpar.begin();
      for(int j = 0; j < i; j++)
          modkey[i]++;
  }
  
  int aij;
    mgood = 0;
    mbad = 0;
    for (int i = 0; i < modules; i++) {
        oc[i] = 0;
        ot[i] = 0;
        whi = *modkey[i];
        if(whi.count(cu_gen)!=0){
            for (ite = whi.begin(); ite != whi.end(); ite++){
                if(*ite >= nodes){
                    cout << "[Error]: Nodes in partition does not belong to network in GraphI::iteractions_to_attractor.\n";
                    exit(1);
                }
                aij = weight(*ite, cu_gen);
                if(aij != 0){
                    if(s[cu_gen]*s[*ite]*aij == 1)
                        mgood++;
                    else
                        mbad++;
                }
            }
                
        }
        else{
            for (ite = whi.begin(); ite != whi.end(); ite++){
                if(*ite >= nodes){
                    cout << "[Error]: Nodes in partition does not belong to network in GraphI::iteractions_to_attractor.\n";
                    exit(1);
                }
                aij = weight(*ite, cu_gen);                                        /*cout << aij << " " << *ite << "\t";*/
                if(aij != 0){
                    if(s[cu_gen]*s[*ite]*aij == 1)
                        oc[i]++;
                    else
                        ot[i]++;
                }
            }
        }
    }
    whi.clear();
    delete[] modkey;
}

void GraphI::set_as_state(int *vec)
{
  int i;
  for (i=0; i< size; i++)
    edo[i] = vec[i];
}

bool GraphI::find_an_attractor_trunc(int steps)
{
  if (yatra) {
    cout << "[Error]: Attractor already found when GraphI::find_an_attractor was called.\n";
    exit(1);
  }
  int i, j, k, l;
  int maxsta = steps;
  int **trajectory;
  trajectory = new int*[maxsta];
  for (i = 0; i < maxsta; i++) 
		trajectory[i] = new int[size];
  k = 0;
  do {
    for (i = 0; i < size; i++)
      trajectory[k][i] = edo[i];
    synchrony();
    k++;
    if (basic.eqvec(edo, size, trajectory[k-1], size)) {
			l = k-1;
			break;
		}
    l = basic.vecinmat(trajectory, k, size, edo, size);
  } while (l==(-1) && k < maxsta);
  if(k >= maxsta){
      for (i = 0; i < maxsta; i++)
        delete [] trajectory[i];
      delete [] trajectory;
      return false;
  }
  else{
    attractor = new int*[k-l];
    for (i= 0; i<(k-l); i++) 
            attractor[i] = new int[size];
    for (i = 0; i < (k-l); i++)
            for (j = 0; j < size; j++)
                attractor[i][j] = trajectory[i+l][j];
    atsize = k-l;
    palen = l;
    for (i = 0; i < maxsta; i++)
        delete [] trajectory[i];
    delete [] trajectory;
    yatra = true;
    return true;
  }
}

void GraphI::synchrony()
{
  int i, j, sum;
  for (i = 0; i < size; i++) {
		sum = 0;
		for (j = 0; j < size; j++) {
			if (weight(j,i) != 0)
        sum = sum + (edo[j]*weight(j,i));
		}
		if (sum > 0)
			ima[i] = 1;
		if (sum < 0)
			ima[i] = -1;
		if (sum == 0)
			ima[i] = edo[i]; //it was 0 before 2017
	}
  for (i = 0; i < size; i++)
    edo[i] = ima[i];
}

int GraphI::attractor_size()
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_size was called.\n";
    exit(1);
  }
  return atsize;
}

int GraphI::attractor_element(int row, int node)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphI::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsize)||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphI::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphI::attractor_element.\n";
  }
  return attractor[row][node];
}

double GraphI::distance_pub(int *goal)
{
    int i,j;
    int dif = 0;
    for (i=0; i < atsize; i++)
        for (j=0; j < size; j++)
        if (goal[j] != attractor[i][j])
            dif++;
    return (1.0 - double(dif)/double(size*atsize));
}

double GraphI::distance_pub(int **goal, int nugs)
{
  
  int numatr = atsize;
  int tam = size; 
  int i,j,k;
  double dis;
  int mulng, mulna;
  if (nugs != numatr) {
    if (nugs < numatr) {
      if ((numatr%nugs) == 0) {
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
          nuatr[(numatr*i)+j][k] = attractor[j][k];
    dis = distance_aux_pub(nugoal, (nugs*mulng), nuatr, (numatr*mulna), tam);
    for (i= 0; i<(nugs*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(numatr*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
   }
  else{
    dis = distance_aux_pub(goal, nugs, attractor, atsize, tam);
  }
  
  return (1.0 - dis);
}


double GraphI::distance_aux_pub(int **goal, int nugs, int **atr, int numatr, int tam)
{
  if (nugs != numatr) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists, d;
  dists = new double[nugs];
  basic.fillv0(dists, nugs);
  for (i=0; i<nugs; i++) {
    dists[i] =0;
    for (j=0; j<nugs; j++) {
      for (k=0;k<tam;k++)
        if (goal[j][k] != atr[(j+i)%nugs][k]){
          dists[i] = dists[i] + (1.0/double(nugs));
        }
    }
  }
  d = basic.find_min(dists, nugs);
  d = d/double(tam);
  delete [] dists;
  return d;
}
//read-print
void GraphI::get_dir_nw_from_file(int nn, string arch)
{
	ifstream sal;
	basic.open_ifstream(sal, arch);
	get_dir_nw_from_file(nn, sal);	
	sal.close();
}

void GraphI::get_dir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, true);
  int i, j, k;
  while (en >> j) {
		en >> i;
		en >> k;
		if (((k==1) || (k==(-1))) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
			nw[i][j] = k;
		else {
			cout << "[Error]: This is not an int network. GraphI::get_dir_nw_from_file.\n";
			exit(1);
		}
	}
  set_default_exclusive_vars();
}

void GraphI::get_dir_nw_from_file_bignw(int nn, istream& en, ifstream& index, int posicion, int elementos)
{
  make_nw(nn, true);
  int i, j, k;
  char ch;
  unsigned long address;
  address = get_address_from_index(index, posicion, elementos);
  
  en.seekg(address, ios::beg);
  while (!en.eof()){
       en >> ch;
       if (ch == '$') break;
  }
  
  while (en >> j && !en.eof()) {
        if (j == 2000) break;
		en >> i;
		en >> k;
		if (((k==1) || (k==(-1))) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
			nw[i][j] = k;
		else {
			cout << "[Error]: This is not an int network. GraphI::get_dir_nw_from_file.\n";
			exit(1);
		}
	}


  if (en.eof()){
			cout << "[Error]: A proihibited address has been selected when trying to acces get_dir_from_file\n";
			exit(1);
		}
  set_default_exclusive_vars();
}

unsigned long GraphI::get_address_from_index(ifstream& inventario, int posicion, int elementos)
{
  unsigned long long lugar;
  string elemento;

  inventario.seekg(0, ios::end);
  unsigned long long count = inventario.tellg();
  
  unsigned long long size = (count)/(unsigned long long)elementos;
  
  inventario.seekg((unsigned long long)posicion*size, ios::beg);
  inventario >> elemento;
  if ("network"+basic.inttostring(posicion) != elemento){
			cout << "[Error]: Could not find the given element in GraphI::get_address_from_index.\n";
			exit(1);
		}
  inventario >> lugar;

  return lugar;  
}

void GraphI::printnw(ostream& sal)
{
  int i, j;
  for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (nw[i][j] == 0)
				sal << 0 << " ";
			else
				if (nw[i][j] == 1)
					sal << "+ ";
				else 
					if (nw[i][j] == -1)
						sal << "- ";
					else {
						cout << "[Error]: error while printing network. GraphI::printnw.\n";
						exit(1);
					}
		}
		sal << endl;
	}
}

void GraphI::export_nw(string arch) 
{
  ofstream fs;
  basic.open_ofstream(fs, arch);
  export_nw(fs);
  fs.close();
}

void GraphI::export_nw(ostream& fs) 
{
  int i, j;
  if (directed) {
    for (i = 0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
  else {
    for (i = 0; i < size; i++)
      for (j=0; j <= i; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
}

void GraphI::export_nw_fixed_width_bignw(ostream& fs, ofstream& key, int element) 
{
  int i, j;
  key << setw(32) << "network"+basic.inttostring(element) << "\t" << setw(32) << fs.tellp() << endl;
  fs << "$" << endl;
  if (directed) {
    for (i = 0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i) != 0)
          fs << setw(6) << j << "\t" << setw(6) << i << "\t" << setw(4) << weight(j, i) << endl;
  }
  else {
    for (i = 0; i < size; i++)
      for (j=0; j <= i; j++)
        if (weight(j,i) != 0)
            fs << setw(6) << j << "\t" << setw(6) << i << "\t" << setw(4) << weight(j, i) << endl;  }
  fs << 2000 << endl << endl << endl;
}

void GraphI::print_dot(string archi)
{
  int i, j;
	ofstream sal;
	string arch;
	arch = archi+".dot";
	basic.open_ofstream(sal, arch);
	string instr, ke;
  if (directed) {
		instr = "digraph G {\n";
		ke = "->";
	}
	else {
		instr = "graph G {\n";
		ke = "--";
	}
	sal << instr;
	if (directed) {
		for (i=0; i < size; i++)
			for (j=0; j < size; j++)
				if (weight(j,i) !=0)
					sal << "\t" << j << "  " << ke << "  " << i << ";\n";
	}
	else 
		for (i = 0; i < size; i++)
			for (j = 0; j <= i; j++)
				if (weight(j,i) !=0)
					sal << "\t" << j << "  " << ke << "  " << i << ";\n";
	sal << "}\n";
	sal.close();
}

void GraphI::print_dot_circ(string archi, double radius)
{
  int i, j;
  ofstream sal;
  string arch;
  arch = archi+".dot";
  basic.open_ofstream(sal, arch);
  string instr, ke,brt,frt;
  if (directed) {
    instr = "digraph G {\n";
    ke = "->";
  }
  else {
    instr = "graph G {\n";
    ke = "--";
  }
  sal << instr;
  sal << "layout = neato;\n";
  sal << "node[shape=circle];\n";
  double radians[size], x[size], y[size], stera, PI=3.141592;
  stera = (2*PI)/(size*1.0);
  radians[0] = 0.75*PI;
  for (i=1; i<size; i++) {
    radians[i] = radians[i-1]+stera;
    if (radians[i] >= (2*PI))
      radians[i] = radians[i] - (2*PI);
  }
  for (i=0; i<size; i++) {
    basic.polar_to_cartesian(radians[i], radius, x[i], y[i]);
    sal << i << "[pos = \"" << x[i] << "," << y[i] << "!\",label=\"" << i << "\"];\n";
  }
  
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++) {
        if ((weight(j,i)) != 0) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  else {
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++) {
        if ((weight(j,i)) != 0) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  sal << "}\n";
  sal.close();
}

void GraphI::print_attractor(ostream& sal)
{
  if (!yatra) {
    cout << "[Error]: Attractor have not been defined when GraphI::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < atsize; i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(i,j) << "\t";
    sal << endl;
  }
}

bool GraphI::robustness_1mutation_gaps(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int nugaps, double &rob)
{
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, nugaps)){
    return false;
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int igual=0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  int nu_mut = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps))
                igual++;
        
            vac.change_interaction(j, i,-1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps))
                igual++;
        }
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps))
                igual += 1;
        }
        
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
      }
  }
  
  vac.clear();
  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
   
  rob = double(igual)/double(nu_mut);
  return true;
}

bool GraphI::robustness_1mutation_cap(int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, double &rob)
{
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s(nodes, mgood, mbad, oc, ot, modules)){
    return false;
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int igual=0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  int nu_mut = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules))
            igual++;
        
            vac.change_interaction(j, i,-1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules))
            igual++;
        }
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules))
            igual += 1;
        }
        
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
      }
  }
  
  vac.clear();
  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
   
  rob = double(igual)/double(nu_mut);
  return true;
}

bool GraphI::rob_by_dist_1mutation_gaps(int**ats, int nugaps, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double &rob, double &dist, double &robprom, double &robgan, double &robper, double &dist_gan, double &dist_per, double &mrob_gan, double &mrob_per, int &count_gan, int &count_per, int &trunc)
{
    
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, nugaps)){
    return false;
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int igual=0;
  double d_at;
  double rpigual = 0;
  int rpigual_ic = 0;
  double digual = 0;
  double igual_ic = 0;
  
  count_gan = 0;
  count_per = 0;
  
  double igual_gan = 0;
  double igual_per = 0;
  double rpigual_gan = 0;
  double rpigual_per = 0;
  double digual_gan = 0;
  double digual_per = 0;
  double digual_ic_div, rpigual_ic_div;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  int count;
  int nu_mut = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                igual++;
                digual += 1;
                rpigual += 1;
                
                igual_gan += 1;
                digual_gan += 1;
                rpigual_gan += 1;
                count_gan++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nugaps; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic++;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                rpigual_gan += (double)rpigual_ic_div/(double)count;
                digual_gan += digual_ic_div/(double)count;
                count_gan++;
            }
        
            vac.change_interaction(j, i,-1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                igual++;
                digual += 1;
                rpigual += 1;
                
                igual_gan += 1;
                digual_gan += 1;
                rpigual_gan += 1;
                count_gan++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nugaps; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic++;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                
                rpigual_gan += (double)rpigual_ic_div/(double)count;
                digual_gan += digual_ic_div/(double)count;
                count_gan++;
            }
        }
        
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                igual = igual + 1;
                digual += 1;
                rpigual += 1;
                
                igual_per += 1;
                digual_per += 1;
                rpigual_per += 1;
                count_per++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nugaps; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic += 1;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                        }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                
                rpigual_per += (double)rpigual_ic_div/(double)count;
                digual_per += digual_ic_div/(double)count;
                count_per++;
            }
        }
        
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
      }
  }
  
  vac.clear();
  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
   
  robprom = double(rpigual)/double(nu_mut);
  rob = double(igual)/double(nu_mut);
  dist = double(digual)/double(nu_mut);
  
  robper = igual_per;
  robgan = igual_gan;
  mrob_gan = rpigual_gan;
  mrob_per = rpigual_per;
  dist_gan = digual_gan;
  dist_per = digual_per;
  return true;
}

bool GraphI::rob_by_dist_1mutation_cap(int**ats, int nic, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double &rob, double &dist, double &robprom,double &robgan, double &robper, double &dist_gan, double &dist_per, double &mrob_gan, double &mrob_per, int &count_gan, int &count_per, int &trunc)
{
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s(nodes, mgood, mbad, oc, ot, modules)){
    return false;
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int igual=0;
  double d_at;
  double rpigual = 0;
  int rpigual_ic = 0;
  double digual = 0;
  double igual_ic = 0;
  
  count_gan = 0;
  count_per = 0;
  
  double igual_gan = 0;
  double igual_per = 0;
  double rpigual_gan = 0;
  double rpigual_per = 0;
  double digual_gan = 0;
  double digual_per = 0;
  double digual_ic_div, rpigual_ic_div;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  int count;
  int nu_mut = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i, 1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                igual++;
                digual += 1;
                rpigual += 1;
                
                igual_gan += 1;
                digual_gan += 1;
                rpigual_gan += 1;
                count_gan++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nic; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic++;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                rpigual_gan += (double)rpigual_ic_div/(double)count;
                digual_gan += digual_ic_div/(double)count;
                count_gan++;
            }
        
            vac.change_interaction(j, i,-1);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                igual++;
                digual += 1;
                rpigual += 1;
                
                igual_gan += 1;
                digual_gan += 1;
                rpigual_gan += 1;
                count_gan++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nic; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic++;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                
                rpigual_gan += (double)rpigual_ic_div/(double)count;
                digual_gan += digual_ic_div/(double)count;
                count_gan++;
            }
        }
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            nu_mut++;
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                igual = igual + 1;
                digual += 1;
                rpigual += 1;
                
                
                igual_per += 1;
                digual_per += 1;
                rpigual_per += 1;
                count_per++;
            }
            else{
                count = 0;
                igual_ic = 0;
                rpigual_ic = 0;
                rpigual_ic_div = 0;
                digual_ic_div = 0;
                for(int ic = 0; ic < nic; ic++){
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at == 1){
                            rpigual_ic += 1;
                            rpigual_ic_div += 1;
                        }
                        igual_ic += d_at;
                        digual_ic_div += d_at;
                        count++;
                        }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
                rpigual += (double)rpigual_ic/(double)count;
                digual += igual_ic/(double)count;
                
                rpigual_per += (double)rpigual_ic_div/(double)count;
                digual_per += digual_ic_div/(double)count;
                count_per++;
            }
        }
        
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, cu_gen, oldval);
      }
  }
  
  vac.clear();
  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
  
  robprom = double(rpigual)/double(nu_mut);   
  rob = double(igual)/double(nu_mut);
  dist = double(digual)/double(nu_mut);
  
  
  robper = igual_per;
  robgan = igual_gan;
  mrob_gan = rpigual_gan;
  mrob_per = rpigual_per;
  dist_gan = digual_gan;
  dist_per = digual_per;
  return true;
}

void GraphI::access_phen_1mutation_gaps(int**ats, int nugaps, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double prop_max_ac, int *tot_pht, double *simpht,  double *modspht, int *pht_g1, int *pht_d1, int *pht_10, int *pht_m1, int *pht_multchan, int *tot_acph, double *simacph, double *modsacph, int *acph_g1, int *acph_d1, int *acph_10, int *acph_m1, int *acph_multchan, int &maxsize, double &phen_size, double &dist_acph, int &trunc)
{
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, nugaps)){
      cout << "[Error]: Network did not produce the modular attractors in GraphI::access_phen_1mutation_gaps.\n";
      exit(1);
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int *nacph_ic;
  nacph_ic = new int[nugaps];
  basic.fillv0(nacph_ic, nugaps);
  int ****ac_phen;
  basic.create_array(ac_phen, nugaps, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph;
  basic.create_array(size_acph, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **repphen;
  basic.create_array(repphen, nugaps, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen, nugaps, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen;
  basic.create_array(sim_phen, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen;
  basic.create_array(mods_phen, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen;
  basic.create_array(totch_phen, nugaps, (int)(prop_max_ac*nodes*nodes));
  
  bool f_nphen;
  
  double *phsize_ic;
  phsize_ic = new double[nugaps];
  basic.fillv0(phsize_ic, nugaps);
  double *dparac;
  dparac = new double[nugaps];
  basic.fillv0(dparac, nugaps);
  double dparac_temp;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  double d_at, d_acph;
  maxsize = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1.0){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
        
            vac.change_interaction(j, i,-1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                 for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
       }
   }
  
  basic.fillv0(simpht, nugaps);
  basic.fillv0(modspht, nugaps);
  basic.fillv0(pht_d1, nugaps);
  basic.fillv0(pht_g1, nugaps);
  basic.fillv0(pht_10, nugaps);
  basic.fillv0(pht_m1, nugaps);
  basic.fillv0(pht_multchan, nugaps);
  
  basic.fillv0(simacph, nugaps);
  basic.fillv0(modsacph, nugaps);
  basic.fillv0(acph_d1, nugaps);
  basic.fillv0(acph_g1, nugaps);
  basic.fillv0(acph_10, nugaps);
  basic.fillv0(acph_m1, nugaps);
  basic.fillv0(acph_multchan, nugaps);
      
  for(int ic = 0; ic < nugaps; ic++){
      tot_pht[ic] = nacph_ic[ic] + basic.sumatoria(repphen[ic], nacph_ic[ic]);
      tot_acph[ic] = nacph_ic[ic];

      for(int i = 0; i < nacph_ic[ic]; i++){
          simpht[ic] += sim_phen[ic][i]*(repphen[ic][i]+1);
          
          if(totch_phen[ic][i] == 1){
              pht_g1[ic] += repphen[ic][i]+1;
              acph_g1[ic]++;
          }
          
          if((1.0-sim_phen[ic][i])*(double)nodes <= 1.00000001){
              pht_d1[ic] += repphen[ic][i]+1;
              acph_d1[ic]++;
          }
          if((1.0-sim_phen[ic][i]) <= 0.1){
              pht_10[ic] += repphen[ic][i]+1;
              acph_10[ic]++;
          }
      }
      if(tot_acph[ic] > 0){
          simacph[ic] = basic.get_mean(sim_phen[ic], tot_acph[ic]);
          simpht[ic] /= double(tot_pht[ic]);
      }
      
      for(int i = 0; i < nacph_ic[ic]; i++){
        if(totch_phen[ic][i] > 1){
            pht_multchan[ic] += repphen[ic][i]+1;
            acph_multchan[ic]++;
            
            modspht[ic] += mods_phen[ic][i]*(repphen[ic][i]+1);
            modsacph[ic] += mods_phen[ic][i];
            
            if(mods_phen[ic][i] == 1){
                pht_m1[ic] += repphen[ic][i]+1;
                acph_m1[ic]++;
            }
        }
      }
      if(acph_multchan[ic] != 0){
          modspht[ic] /= double(pht_multchan[ic]);
          modsacph[ic] /= double(acph_multchan[ic]);
      }
      
      
      if(tot_pht[ic]!= 0)
        phsize_ic[ic] = phsize_ic[ic]/double(tot_pht[ic]);

      if(tot_acph[ic] > 1)
        dparac[ic] = dparac[ic]/(double(tot_acph[ic]*(tot_acph[ic] - 1))/2.0);
  }
  
  dist_acph = basic.get_mean(dparac, nugaps);
  phen_size = basic.get_mean(phsize_ic, nugaps);
  
  
  vac.clear();
  for(int i = 0; i < nugaps; i++){
    for(int j = 0; j < (int)(prop_max_ac*nodes*nodes); j++){
        for(int k = 0; k < at_step; k++){
            delete ac_phen[i][j][k];
        }
        delete ac_phen[i][j];
    }
    delete ac_phen[i];
    delete size_acph[i];
    delete repphen[i];
    delete sim_phen[i];
    delete mods_phen[i];
    delete totch_phen[i];
  }
  delete[] nacph_ic;
  delete[] ac_phen;
  delete[] size_acph;
  delete[] repphen;
  delete[] sim_phen;
  delete[] mods_phen;
  delete[] totch_phen;
  delete[] phsize_ic;
  delete[] dparac;

  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
}

void GraphI::access_phen_1mutation_cap(int**ats, int nuic, int *s, int* mgood, int* mbad, int** oc, int** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double prop_max_ac, int *tot_pht, double *simpht,  double *modspht, int *pht_g1, int *pht_d1, int *pht_m1, int *pht_multchan, int *tot_acph, double *simacph, double *modsacph, int *acph_g1, int *acph_d1, int *acph_m1, int *acph_multchan, int &maxsize, double &phen_size, double &dist_acph, int &trunc)
{
  interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
  FitnessI law(est);
  
  if(!law.cap_to_s(nodes, mgood, mbad, oc, ot, modules)){
      cout << "[Error]: Network did not produce the modular attractors in GraphI::access_phen_1mutation_gaps.\n";
      exit(1);
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[i];
      c_mbad[i] = mbad[i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[i][j];
        c_ot[i][j] = ot[i][j];
      }
  }
  
  int *nacph_ic;
  nacph_ic = new int[nuic];
  basic.fillv0(nacph_ic, nuic);
  int ****ac_phen;
  basic.create_array(ac_phen, nuic, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph;
  basic.create_array(size_acph, nuic, (int)(prop_max_ac*nodes*nodes));
  int **repphen;
  basic.create_array(repphen, nuic, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen, nuic, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen;
  basic.create_array(sim_phen, nuic, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen;
  basic.create_array(mods_phen, nuic, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen;
  basic.create_array(totch_phen, nuic, (int)(prop_max_ac*nodes*nodes));
  
  bool f_nphen;
  
  double *phsize_ic;
  phsize_ic = new double[nuic];
  basic.fillv0(phsize_ic, nuic);
  double *dparac;
  dparac = new double[nuic];
  basic.fillv0(dparac, nuic);
  double dparac_temp;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  double d_at, d_acph;
  maxsize = 0;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1.0){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
        
            vac.change_interaction(j, i,-1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                 for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1.0){
                            phsize_ic[ic] += vac.attractor_size();
                            if(vac.attractor_size() > maxsize)
                                maxsize = vac.attractor_size();
                            dparac_temp = 0;
                            for(int nac = 0; nac < nacph_ic[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen[ic][nac], size_acph[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen[ic][nac]++;
                                    break;
                                 }
                                 dparac_temp += d_acph;
                             }
                             if(f_nphen){
                                mods_phen[ic][nacph_ic[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen[ic][nacph_ic[ic]], sim_phen[ic][nacph_ic[ic]]);
                                size_acph[ic][nacph_ic[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph[ic][nacph_ic[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen[ic][nacph_ic[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                dparac[ic] += dparac_temp;
                                nacph_ic[ic]++;
                                if(nacph_ic[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::access_phen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[cu_gen];
        c_mbad[cu_gen] = mbad[cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
       }
   }
  
  basic.fillv0(simpht, nuic);
  basic.fillv0(modspht, nuic);
  basic.fillv0(pht_g1, nuic);
  basic.fillv0(pht_d1, nuic);
  basic.fillv0(pht_m1, nuic);
  basic.fillv0(pht_multchan, nuic);
  
  basic.fillv0(simacph, nuic);
  basic.fillv0(modsacph, nuic);
  basic.fillv0(acph_g1, nuic);
  basic.fillv0(acph_d1, nuic);
  basic.fillv0(acph_m1, nuic);
  basic.fillv0(acph_multchan, nuic);
      
  for(int ic = 0; ic < nuic; ic++){
      tot_pht[ic] = nacph_ic[ic] + basic.sumatoria(repphen[ic], nacph_ic[ic]);
      tot_acph[ic] = nacph_ic[ic];

      for(int i = 0; i < nacph_ic[ic]; i++){
          simpht[ic] += sim_phen[ic][i]*(repphen[ic][i]+1);
          
          if(totch_phen[ic][i] == 1){
              pht_g1[ic] += repphen[ic][i]+1;
              acph_g1[ic]++;
          }
          
          if((1.0-sim_phen[ic][i])*(double)nodes <= 1.00000001){
              pht_d1[ic] += repphen[ic][i]+1;
              acph_d1[ic]++;
          }
      }
      if(nacph_ic[ic] > 0){
          simacph[ic] = basic.get_mean(sim_phen[ic], tot_acph[ic]);
          simpht[ic] /= tot_pht[ic];
      }
      
      for(int i = 0; i < nacph_ic[ic]; i++){
        if(totch_phen[ic][i] > 1){
            pht_multchan[ic] += repphen[ic][i]+1;
            acph_multchan[ic]++;
            
            modspht[ic] += mods_phen[ic][i]*(repphen[ic][i]+1);
            modsacph[ic] += mods_phen[ic][i];
            
            if(mods_phen[ic][i] == 1){
                pht_m1[ic] += repphen[ic][i]+1;
                acph_m1[ic]++;
            }
        }
      }
      if(acph_multchan[ic] != 0){
          modspht[ic] /= double(pht_multchan[ic]);
          modsacph[ic] /= double(acph_multchan[ic]);
      }
      
      
      if(tot_pht[ic]!= 0)
        phsize_ic[ic] = phsize_ic[ic]/double(tot_pht[ic]);

      if(tot_acph[ic] > 1)
        dparac[ic] = dparac[ic]/((double)(tot_acph[ic]*(tot_acph[ic] - 1))/2.0);
  }
  
  dist_acph = basic.get_mean(dparac, nuic);
  phen_size = basic.get_mean(phsize_ic, nuic);
  
  
  vac.clear();
  for(int i = 0; i < nuic; i++){
    for(int j = 0; j < (int)(prop_max_ac*nodes*nodes); j++){
        for(int k = 0; k < at_step; k++){
            delete ac_phen[i][j][k];
        }
        delete ac_phen[i][j];
    }
    delete ac_phen[i];
    delete size_acph[i];
    delete repphen[i];
    delete sim_phen[i];
    delete mods_phen[i];
    delete totch_phen[i];
  }
  delete[] nacph_ic;
  delete[] ac_phen;
  delete[] size_acph;
  delete[] repphen;
  delete[] sim_phen;
  delete[] mods_phen;
  delete[] totch_phen;
  delete[] phsize_ic;
  delete[] dparac;

  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
}

void GraphI::access_phen_translape_gaps(GraphI &temp, int**ats, int nugaps, int *s, int** mgood, int** mbad, int*** oc, int*** ot, int nodes, int modules, set<set<int> >& predpar, int ***intmods, int at_step, double prop_max_ac, int **tot_phen, int **tot_acphen, int *tot_acphens, double *sim_acphens,  double *mods_acphens, int *acphens_d1, int *acphens_m1, int *acphens_multchan, int *tot_newpht, double *simnewpht,  double *modsnewpht, int *newpht_g1, int *newpht_d1, int *newpht_m1, int *newpht_multchan, int *tot_newacph, double *simnewacph, double *modsnewacph, int *newacph_g1, int *newacph_d1, int *newacph_m1, int *newacph_multchan, int &trunc)
{
  if(temp.size != size){
    cout << "[Error]: nw have different size in GraphI::access_phen_translape_cap.\n";
    exit(1);
  }
  FitnessI law(est);
  
  interactions_to_attractor(s, mgood[0], mbad[0], oc[0], ot[0], nodes, modules, predpar);
  if(!law.cap_to_s_gaps(nodes, mgood[0], mbad[0], oc[0], ot[0], modules, intmods, nugaps)){
      cout << "[Error]: First network did not produce the modular attractors in GraphI::access_phen_translape_gaps.\n";
      exit(1);
  }
  temp.interactions_to_attractor(s, mgood[1], mbad[1], oc[1], ot[1], nodes, modules, predpar);
  if(!law.cap_to_s_gaps(nodes, mgood[1], mbad[1], oc[1], ot[1], modules, intmods, nugaps)){
      cout << "[Error]: Second network did not produce the modular attractors in GraphI::access_phen_translape_gaps.\n";
      exit(1);
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[0][i];
      c_mbad[i] = mbad[0][i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[0][i][j];
        c_ot[i][j] = ot[0][i][j];
      }
  }
  
  int *nacph_ic1;
  nacph_ic1 = new int[nugaps];
  basic.fillv0(nacph_ic1, nugaps);
  int ****ac_phen1;
  basic.create_array(ac_phen1, nugaps, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph1;
  basic.create_array(size_acph1, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **repphen1;
  basic.create_array(repphen1, nugaps, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen1, nugaps, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen1;
  basic.create_array(sim_phen1, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen1;
  basic.create_array(mods_phen1, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen1;
  basic.create_array(totch_phen1, nugaps, (int)(prop_max_ac*nodes*nodes));
  
  bool f_nphen;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  double d_at, d_acph;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
        
            vac.change_interaction(j, i,-1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                 for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[0][cu_gen];
        c_mbad[cu_gen] = mbad[0][cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[0][cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[0][cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
       }
   }
   
   
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[1][i];
      c_mbad[i] = mbad[1][i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[1][i][j];
        c_ot[i][j] = ot[1][i][j];
      }
  }
  
  GraphI temp_vac(est);
  int *nacph_ic2;
  nacph_ic2 = new int[nugaps];
  basic.fillv0(nacph_ic2, nugaps);
  int ****ac_phen2;
  basic.create_array(ac_phen2, nugaps, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph2;
  basic.create_array(size_acph2, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **repphen2;
  basic.create_array(repphen2, nugaps, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen2, nugaps, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen2;
  basic.create_array(sim_phen2, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen2;
  basic.create_array(mods_phen2, nugaps, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen2;
  basic.create_array(totch_phen2, nugaps, (int)(prop_max_ac*nodes*nodes));
  temp.make_copy(temp_vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = temp_vac.weight(j, i);
        if (oldval == 0){
            temp_vac.change_interaction(j, i,1);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
            }
        
            temp_vac.change_interaction(j, i,-1);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            temp_vac.change_interaction(j, i, 0);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen_gaps(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules, intmods[cu_gen], nugaps)){
                 for(int ic = 0; ic < nugaps; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[1][cu_gen];
        c_mbad[cu_gen] = mbad[1][cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[1][cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[1][cu_gen][cumod];
        }
        temp_vac.change_interaction(j, i, oldval);;
       }
  }
   
  basic.fillv0(tot_phen[0], nugaps);
  basic.fillv0(tot_phen[1], nugaps); 
  
  for(int ic = 0; ic < nugaps; ic++){
    for(int i = 0; i < nacph_ic1[ic]; i++){
        tot_phen[0][ic] += repphen1[ic][i] + 1;
    }
    for(int j = 0; j < nacph_ic2[ic]; j++){
        tot_phen[1][ic] += repphen2[ic][j] + 1;
    }
  }
  
  basic.fillv0(tot_acphen[0], nugaps);
  basic.fillv0(tot_acphen[1], nugaps);
    
  bool **phen2in1;
  int maxph2 = basic.find_max(nacph_ic2, nugaps);
  basic.create_array(phen2in1, nugaps, maxph2);
  basic.fillmat0(phen2in1, nugaps, maxph2);
  double dist_phen;
  
  bool f_trans;
  for(int ic = 0; ic < nugaps; ic++){
    for(int i = 0; i < nacph_ic1[ic]; i++){
        f_trans = false;
        for(int j = 0; j < nacph_ic2[ic]; j++){
            dist_phen = law.distance(ac_phen1[ic][i], size_acph1[ic][i], ac_phen2[ic][j], size_acph2[ic][j], size);
            if(dist_phen == 0){
                if(f_trans){
                    cout << "Error two equal phenotypes in accumulated phenotypes in  GraphI::access_phen_translape_cap.\n";
                    exit(1);
                }
                phen2in1[ic][j] = true;
                f_trans = true;
            }        
        }
    }
    tot_acphen[0][ic] = nacph_ic1[ic];
    tot_acphen[1][ic] = nacph_ic2[ic];
  }
  
  basic.fillv0(tot_acphens, nugaps);
  basic.fillv0(sim_acphens, nugaps);
  basic.fillv0(mods_acphens, nugaps);
  basic.fillv0(acphens_d1, nugaps);
  basic.fillv0(acphens_m1, nugaps);
  basic.fillv0(acphens_multchan, nugaps);
  
  for(int ic = 0; ic < nugaps; ic++){

      for(int i = 0; i < nacph_ic1[ic]; i++){
          tot_acphens[ic]++; 
          sim_acphens[ic] += sim_phen1[ic][i];
          if((1.0-sim_phen1[ic][i])*(double)nodes <= 1.00000001)
              acphens_d1[ic]++;
          
          if(totch_phen1[ic][i] > 1){
              acphens_multchan[ic]++;
              mods_acphens[ic] += mods_phen1[ic][i];
              if(mods_phen1[ic][i] == 1)
                  acphens_m1[ic]++;
          }
      }
      
      for(int j = 0; j < nacph_ic2[ic]; j++){
          if(!phen2in1[ic][j]){
              tot_acphens[ic]++; 
              sim_acphens[ic] += sim_phen2[ic][j];
              if((1.0-sim_phen2[ic][j])*(double)nodes <= 1.00000001)
                  acphens_d1[ic]++;
              
              if(totch_phen2[ic][j] > 1){
                  acphens_multchan[ic]++;
                  mods_acphens[ic] += mods_phen2[ic][j];
                  if(mods_phen2[ic][j] == 1)
                      acphens_m1[ic]++;
              }
          }
      }
      if(tot_acphens[ic] > 0)
          sim_acphens[ic] /= double(tot_acphens[ic]);

      if(acphens_multchan[ic] > 0)
          mods_acphens[ic] /= double(acphens_multchan[ic]);
  }
  
  basic.fillv0(tot_newpht, nugaps);
  basic.fillv0(simnewpht, nugaps);
  basic.fillv0(modsnewpht, nugaps);
  basic.fillv0(newpht_g1, nugaps);
  basic.fillv0(newpht_d1, nugaps);
  basic.fillv0(newpht_m1, nugaps);
  basic.fillv0(newpht_multchan, nugaps);
  
  basic.fillv0(tot_newacph, nugaps);
  basic.fillv0(simnewacph, nugaps);
  basic.fillv0(modsnewacph, nugaps);
  basic.fillv0(newacph_g1, nugaps);
  basic.fillv0(newacph_d1, nugaps);
  basic.fillv0(newacph_m1, nugaps);
  basic.fillv0(newacph_multchan, nugaps);
  
  for(int ic = 0; ic < nugaps; ic++){
    for(int j = 0; j < nacph_ic2[ic]; j++){
        if(!phen2in1[ic][j]){
            tot_newpht[ic] += repphen2[ic][j]+1;
            tot_newacph[ic]++;
            
            simnewpht[ic] += sim_phen2[ic][j]*(repphen2[ic][j]+1);
            simnewacph[ic] += sim_phen2[ic][j];
            
            if(totch_phen2[ic][j] == 1){
                newpht_g1[ic] += repphen2[ic][j]+1;
                newacph_g1[ic]++;
            }
                
            if((1.0-sim_phen2[ic][j])*(double)nodes <= 1.00000001){
                newpht_d1[ic] += repphen2[ic][j]+1;
                newacph_d1[ic]++;
            }
          
            if(totch_phen2[ic][j] > 1){
                newpht_multchan[ic] += repphen2[ic][j]+1;
                newacph_multchan[ic]++;
            
                modsnewpht[ic] += mods_phen2[ic][j]*(repphen2[ic][j]+1);
                modsnewacph[ic] += mods_phen2[ic][j];
            
                if(mods_phen2[ic][j] == 1){
                    newpht_m1[ic] += repphen2[ic][j]+1;
                    newacph_m1[ic]++;
                }
            }
        }
    
    }
    if(tot_newacph[ic] > 0){
        simnewacph[ic] /= double(tot_newacph[ic]);
        simnewpht[ic] /= double(tot_newpht[ic]);
    }
    
    if(newacph_multchan[ic] > 0){
        modsnewpht[ic] /= double(newpht_multchan[ic]);
        modsnewacph[ic] /= double(newacph_multchan[ic]);
    }
    
  }
  
  
  vac.clear();
  temp_vac.clear();
  for(int i = 0; i < nugaps; i++){
    for(int j = 0; j < (int)(prop_max_ac*nodes*nodes); j++){
        for(int k = 0; k < at_step; k++){
            delete ac_phen1[i][j][k];
            delete ac_phen2[i][j][k];
        }
        delete ac_phen1[i][j];
        delete ac_phen2[i][j];
    }
    delete ac_phen1[i];
    delete size_acph1[i];
    delete repphen1[i];
    delete sim_phen1[i];
    delete mods_phen1[i];
    delete totch_phen1[i];
    
    delete ac_phen2[i];
    delete size_acph2[i];
    delete repphen2[i];
    delete sim_phen2[i];
    delete mods_phen2[i];
    delete totch_phen2[i];
    
    delete phen2in1[i];
  }
  delete[] nacph_ic1;
  delete[] ac_phen1;
  delete[] size_acph1;
  delete[] repphen1;
  delete[] sim_phen1;
  delete[] mods_phen1;
  delete[] totch_phen1;
    
  delete[] nacph_ic2;
  delete[] ac_phen2;
  delete[] size_acph2;
  delete[] repphen2;
  delete[] sim_phen2;
  delete[] mods_phen2;
  delete[] totch_phen2;
  
  delete[] phen2in1;
  

  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
  
}

void GraphI::access_phen_translape_cap(GraphI &temp, int**ats, int nuic, int *s, int** mgood, int** mbad, int*** oc, int*** ot, int nodes, int modules, set<set<int> >& predpar, int at_step, double prop_max_ac, int **tot_phen, int **tot_acphen, int *tot_acphens, double *sim_acphens,  double *mods_acphens, int *acphens_d1, int *acphens_m1, int *acphens_multchan, int *tot_newpht, double *simnewpht,  double *modsnewpht, int *newpht_g1, int *newpht_d1, int *newpht_m1, int *newpht_multchan, int *tot_newacph, double *simnewacph, double *modsnewacph, int *newacph_g1, int *newacph_d1, int *newacph_m1, int *newacph_multchan, int &trunc)
{
  if(temp.size != size){
    cout << "[Error]: nw have different size in GraphI::access_phen_translape_cap.\n";
    exit(1);
  }
  FitnessI law(est);
  
  interactions_to_attractor(s, mgood[0], mbad[0], oc[0], ot[0], nodes, modules, predpar);
  if(!law.cap_to_s(nodes, mgood[0], mbad[0], oc[0], ot[0], modules)){
      cout << "[Error]: First network did not produce the modular attractors in GraphI::access_phen_translape_gaps.\n";
      exit(1);
  }
  temp.interactions_to_attractor(s, mgood[1], mbad[1], oc[1], ot[1], nodes, modules, predpar);
  if(!law.cap_to_s(nodes, mgood[1], mbad[1], oc[1], ot[1], modules)){
      cout << "[Error]: Second network did not produce the modular attractors in GraphI::access_phen_translape_gaps.\n";
      exit(1);
  }
  
  int *c_mgood, *c_mbad;
  int **c_oc, **c_ot;
  
  c_mgood = new int[nodes];
  c_mbad = new int[nodes];
  c_oc = new int*[nodes];
  c_ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      c_oc[i] = new int[modules];
      c_ot[i] = new int[modules];
  }
  
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[0][i];
      c_mbad[i] = mbad[0][i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[0][i][j];
        c_ot[i][j] = ot[0][i][j];
      }
  }
  
  int *nacph_ic1;
  nacph_ic1 = new int[nuic];
  basic.fillv0(nacph_ic1, nuic);
  int ****ac_phen1;
  basic.create_array(ac_phen1, nuic, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph1;
  basic.create_array(size_acph1, nuic, (int)(prop_max_ac*nodes*nodes));
  int **repphen1;
  basic.create_array(repphen1, nuic, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen1, nuic, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen1;
  basic.create_array(sim_phen1, nuic, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen1;
  basic.create_array(mods_phen1, nuic, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen1;
  basic.create_array(totch_phen1, nuic, (int)(prop_max_ac*nodes*nodes));
  
  bool f_nphen;
  
  trunc = 0;
  GraphI vac(est);
  int cu_gen;
  int oldval;
  double d_at, d_acph;
  make_copy(vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = vac.weight(j, i);
        if (oldval == 0){
            vac.change_interaction(j, i,1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at < 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
        
            vac.change_interaction(j, i,-1);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            vac.change_interaction(j, i, 0);
            
            vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                 for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    vac.set_as_state(ats[ic]);
                    if(vac.find_an_attractor_trunc(at_step)){
                        d_at = vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic1[ic]; nac++){
                                d_acph = vac.distance_pub(ac_phen1[ic][nac], size_acph1[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen1[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen1[ic][nacph_ic1[ic]] = vac.how_many_perturbed_modules(ats[ic], modules, totch_phen1[ic][nacph_ic1[ic]], sim_phen1[ic][nacph_ic1[ic]]);
                                size_acph1[ic][nacph_ic1[ic]] = vac.attractor_size();
                                for(int l = 0; l < size_acph1[ic][nacph_ic1[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen1[ic][nacph_ic1[ic]][l][m] = vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic1[ic]++;
                                if(nacph_ic1[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[0][cu_gen];
        c_mbad[cu_gen] = mbad[0][cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[0][cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[0][cu_gen][cumod];
        }
        vac.change_interaction(j, i, oldval);
       }
   }
   
   
  for(int i = 0; i < nodes; i++){
      c_mgood[i] = mgood[1][i];
      c_mbad[i] = mbad[1][i];
      for(int j = 0; j < modules; j++){
        c_oc[i][j] = oc[1][i][j];
        c_ot[i][j] = ot[1][i][j];
      }
  }
  
  GraphI temp_vac(est);
  int *nacph_ic2;
  nacph_ic2 = new int[nuic];
  basic.fillv0(nacph_ic2, nuic);
  int ****ac_phen2;
  basic.create_array(ac_phen2, nuic, (int)(prop_max_ac*nodes*nodes), at_step, nodes);
  int **size_acph2;
  basic.create_array(size_acph2, nuic, (int)(prop_max_ac*nodes*nodes));
  int **repphen2;
  basic.create_array(repphen2, nuic, (int)(prop_max_ac*nodes*nodes));
  basic.fillmat0(repphen2, nuic, (int)(prop_max_ac*nodes*nodes));
  
  double **sim_phen2;
  basic.create_array(sim_phen2, nuic, (int)(prop_max_ac*nodes*nodes));
  int **mods_phen2;
  basic.create_array(mods_phen2, nuic, (int)(prop_max_ac*nodes*nodes));
  int **totch_phen2;
  basic.create_array(totch_phen2, nuic, (int)(prop_max_ac*nodes*nodes));
  temp.make_copy(temp_vac);
  
  for (int i=0; i<nodes;i++){
      cu_gen = i;
      for (int j=0; j<nodes; j++){
        
        oldval = temp_vac.weight(j, i);
        if (oldval == 0){
            temp_vac.change_interaction(j, i,1);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
            }
        
            temp_vac.change_interaction(j, i,-1);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
          }
        }
         
        else if(oldval != 0){
            temp_vac.change_interaction(j, i, 0);
            
            temp_vac.interactions_to_attractor_1gen(s, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], nodes, modules, predpar, cu_gen);
        
            if(!law.cap_to_s_1gen(cu_gen, c_mgood[cu_gen], c_mbad[cu_gen], c_oc[cu_gen], c_ot[cu_gen], modules)){
                 for(int ic = 0; ic < nuic; ic++){
                    f_nphen = true;
                    temp_vac.set_as_state(ats[ic]);
                    if(temp_vac.find_an_attractor_trunc(at_step)){
                        d_at = temp_vac.distance_pub(ats[ic]);
                        if(d_at != 1){
                            for(int nac = 0; nac < nacph_ic2[ic]; nac++){
                                d_acph = temp_vac.distance_pub(ac_phen2[ic][nac], size_acph2[ic][nac]);
                                if(d_acph == 1){
                                    f_nphen = false;
                                    repphen2[ic][nac]++;
                                    break;
                                 }
                             }
                             if(f_nphen){
                                mods_phen2[ic][nacph_ic2[ic]] = temp_vac.how_many_perturbed_modules(ats[ic], modules, totch_phen2[ic][nacph_ic2[ic]], sim_phen2[ic][nacph_ic2[ic]]);
                                size_acph2[ic][nacph_ic2[ic]] = temp_vac.attractor_size();
                                for(int l = 0; l < size_acph2[ic][nacph_ic2[ic]]; l++){
                                    for(int m = 0; m < size; m++){
                                        ac_phen2[ic][nacph_ic2[ic]][l][m] = temp_vac.attractor_element(l,m);
                                    }
                                }
                                nacph_ic2[ic]++;
                                if(nacph_ic2[ic] >= (int)(prop_max_ac*nodes*nodes)){
                                    cout << "[Error]: Accumulated phenotypes exceed the number of proposed phenotypes in prop_max_ac in GraphI::accessphen_1mutation_cap.\n";
		                            exit(1);
                                }
                            }
                        }
                    }
                    else
                        trunc++;
                    temp_vac.clear_attractor();
                }
            }
         }
         
        c_mgood[cu_gen] = mgood[1][cu_gen];
        c_mbad[cu_gen] = mbad[1][cu_gen];
        for(int cumod = 0; cumod < modules; cumod++){
            c_oc[cu_gen][cumod] = oc[1][cu_gen][cumod];
            c_ot[cu_gen][cumod] = ot[1][cu_gen][cumod];
        }
        temp_vac.change_interaction(j, i, oldval);;
       }
  }
   
  basic.fillv0(tot_phen[0], nuic);
  basic.fillv0(tot_phen[1], nuic); 
  
  for(int ic = 0; ic < nuic; ic++){
    for(int i = 0; i < nacph_ic1[ic]; i++){
        tot_phen[0][ic] += repphen1[ic][i] + 1;
    }
    for(int j = 0; j < nacph_ic2[ic]; j++){
        tot_phen[1][ic] += repphen2[ic][j] + 1;
    }
  }
  
  basic.fillv0(tot_acphen[0], nuic);
  basic.fillv0(tot_acphen[1], nuic);
    
  bool **phen2in1;
  int maxph2 = basic.find_max(nacph_ic2, nuic);
  basic.create_array(phen2in1, nuic, maxph2);
  basic.fillmat0(phen2in1, nuic, maxph2);
  double dist_phen;
  
  bool f_trans;
  for(int ic = 0; ic < nuic; ic++){
    for(int i = 0; i < nacph_ic1[ic]; i++){
        f_trans = false;
        for(int j = 0; j < nacph_ic2[ic]; j++){
            dist_phen = law.distance(ac_phen1[ic][i], size_acph1[ic][i], ac_phen2[ic][j], size_acph2[ic][j], size);
            if(dist_phen == 0){
                if(f_trans){
                    cout << "Error two equal phenotypes in accumulated phenotypes in  GraphI::access_phen_translape_cap.\n";
                    exit(1);
                }
                phen2in1[ic][j] = true;
                f_trans = true;
            }        
        }
    }
    tot_acphen[0][ic] = nacph_ic1[ic];
    tot_acphen[1][ic] = nacph_ic2[ic];
  }
  
  basic.fillv0(tot_acphens, nuic);
  basic.fillv0(sim_acphens, nuic);
  basic.fillv0(mods_acphens, nuic);
  basic.fillv0(acphens_d1, nuic);
  basic.fillv0(acphens_m1, nuic);
  basic.fillv0(acphens_multchan, nuic);
  
  for(int ic = 0; ic < nuic; ic++){

      for(int i = 0; i < nacph_ic1[ic]; i++){
          tot_acphens[ic]++; 
          sim_acphens[ic] += sim_phen1[ic][i];
          if((1.0-sim_phen1[ic][i])*(double)nodes <= 1.00000001)
              acphens_d1[ic]++;
          
          if(totch_phen1[ic][i] > 1){
              acphens_multchan[ic]++;
              mods_acphens[ic] += mods_phen1[ic][i];
              if(mods_phen1[ic][i] == 1)
                  acphens_m1[ic]++;
          }
      }
      
      for(int j = 0; j < nacph_ic2[ic]; j++){
          if(!phen2in1[ic][j]){
              tot_acphens[ic]++; 
              sim_acphens[ic] += sim_phen2[ic][j];
              if((1.0-sim_phen2[ic][j])*(double)nodes <= 1.00000001)
                  acphens_d1[ic]++;
              
              if(totch_phen2[ic][j] > 1){
                  acphens_multchan[ic]++;
                  mods_acphens[ic] += mods_phen2[ic][j];
                  if(mods_phen2[ic][j] == 1)
                      acphens_m1[ic]++;
              }
          }
      }
      if(tot_acphens[ic] > 0)
          sim_acphens[ic] /= double(tot_acphens[ic]);

      if(acphens_multchan[ic] > 0)
          mods_acphens[ic] /= double(acphens_multchan[ic]);
  }
  
  basic.fillv0(tot_newpht, nuic);
  basic.fillv0(simnewpht, nuic);
  basic.fillv0(modsnewpht, nuic);
  basic.fillv0(newpht_g1, nuic);
  basic.fillv0(newpht_d1, nuic);
  basic.fillv0(newpht_m1, nuic);
  basic.fillv0(newpht_multchan, nuic);
  
  basic.fillv0(tot_newacph, nuic);
  basic.fillv0(simnewacph, nuic);
  basic.fillv0(modsnewacph, nuic);
  basic.fillv0(newacph_g1, nuic);
  basic.fillv0(newacph_d1, nuic);
  basic.fillv0(newacph_m1, nuic);
  basic.fillv0(newacph_multchan, nuic);
  
  for(int ic = 0; ic < nuic; ic++){
    for(int j = 0; j < nacph_ic2[ic]; j++){
        if(!phen2in1[ic][j]){
            tot_newpht[ic] += repphen2[ic][j]+1;
            tot_newacph[ic]++;
            
            simnewpht[ic] += sim_phen2[ic][j]*(repphen2[ic][j]+1);
            simnewacph[ic] += sim_phen2[ic][j];
            
            if(totch_phen2[ic][j] == 1){
                newpht_g1[ic] += repphen2[ic][j]+1;
                newacph_g1[ic]++;
            }
            
            if((1.0-sim_phen2[ic][j])*(double)nodes <= 1.00000001){
                newpht_d1[ic] += repphen2[ic][j]+1;
                newacph_d1[ic]++;
            }
          
            if(totch_phen2[ic][j] > 1){
                newpht_multchan[ic] += repphen2[ic][j]+1;
                newacph_multchan[ic]++;
            
                modsnewpht[ic] += mods_phen2[ic][j]*(repphen2[ic][j]+1);
                modsnewacph[ic] += mods_phen2[ic][j];
            
                if(mods_phen2[ic][j] == 1){
                    newpht_m1[ic] += repphen2[ic][j]+1;
                    newacph_m1[ic]++;
                }
            }
        }
    
    }
    if(tot_newacph[ic] > 0){
        simnewacph[ic] /= double(tot_newacph[ic]);
        simnewpht[ic] /= double(tot_newpht[ic]);
    }
    
    if(newacph_multchan[ic] > 0){
        modsnewpht[ic] /= double(newpht_multchan[ic]);
        modsnewacph[ic] /= double(newacph_multchan[ic]);
    }
    
  }
  
  
  vac.clear();
  temp_vac.clear();
  for(int i = 0; i < nuic; i++){
    for(int j = 0; j < (int)(prop_max_ac*nodes*nodes); j++){
        for(int k = 0; k < at_step; k++){
            delete ac_phen1[i][j][k];
            delete ac_phen2[i][j][k];
        }
        delete ac_phen1[i][j];
        delete ac_phen2[i][j];
    }
    delete ac_phen1[i];
    delete size_acph1[i];
    delete repphen1[i];
    delete sim_phen1[i];
    delete mods_phen1[i];
    delete totch_phen1[i];
    
    delete ac_phen2[i];
    delete size_acph2[i];
    delete repphen2[i];
    delete sim_phen2[i];
    delete mods_phen2[i];
    delete totch_phen2[i];
    
    delete phen2in1[i];
  }
  delete[] nacph_ic1;
  delete[] ac_phen1;
  delete[] size_acph1;
  delete[] repphen1;
  delete[] sim_phen1;
  delete[] mods_phen1;
  delete[] totch_phen1;
    
  delete[] nacph_ic2;
  delete[] ac_phen2;
  delete[] size_acph2;
  delete[] repphen2;
  delete[] sim_phen2;
  delete[] mods_phen2;
  delete[] totch_phen2;
  
  delete[] phen2in1;
  

  for(int i = 0; i < nodes; i++){
      delete c_oc[i];
      delete c_ot[i];
  }
  delete[] c_mgood;
  delete[] c_mbad;
  delete[] c_oc;
  delete[] c_ot;
  
}

int GraphI::how_many_perturbed_modules(int *attorig, int modules, int& total_chan, double &sim)
{
  total_chan = 0;
  int dif_mod = 0;
  int attam;
  int **unat;
  int tam = number_of_nodes();
  if (tam%modules != 0){
        cout << "[Error]: nodes in network are not a multiple integer of modules in GraphI::which_module \n";
        exit(1);
    }
  int sep = tam/modules;
  bool ya;
  
    attam = 0;  
    attam = attractor_size();
    
    basic.create_array(unat, attam, tam);
    for (int j = 0; j < attam; j++)
        for (int k = 0; k < tam; k++)
            unat[j][k] = attractor_element(j, k);
    
    for(int cu_mod = 0; cu_mod < modules; cu_mod++){
        ya = false;
        for(int j = 0; j < attam; j++){
            for(int n = 0; n < sep; n++){
                if(attorig[sep*cu_mod + n] != unat[j][sep*cu_mod + n]){
                    total_chan++;
                    if(!ya){
                        dif_mod++;
                        ya = true;
                    }
                }
            }
        }
    }
    
    sim = 1.0 - double(total_chan)/double(attam*tam);
    
    for (int j = 0; j < attam; j++)
        delete unat[j];
    delete[] unat;
    
    return dif_mod;
}
//private

void GraphI::set_default_exclusive_vars()
{
  yae = false;
  yadeg = false;
  yadegdist=false;
  yatra = false;
  yatras = false;
  yamadya = false;
}

void GraphI::prepare_for_degrees()
{
	degree = new int[size];
	outdegree = new int[size];
	indegree = new int[size];
}

void GraphI::prepare_for_degdist()
{
	degdist = new int[size+1];
	basic.fillv0(degdist, size+1);
	
	odegdist = new int[size+1];
	basic.fillv0(odegdist, size+1);
	
	idegdist = new int[size+1];
	basic.fillv0(idegdist, size+1);
}
