#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
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
using namespace std;

int main()
{

  ofstream fs;
  ofstream omuestra;
  ofstream oindex;
  int seed = 24125;
  int nodes = 24;
  int modules = 4;
  int samplesize = 5000;
  int stepsize = 20*nodes*nodes;
  double conectivity = (double)3/(double)nodes;
  int genpermod = nodes/modules;
  int nugaps = modules;
  
  int wol = (int)nodes*nodes*conectivity;
  Alea jacta(seed);
  Basics bas(jacta);
  Cap_evolI cap(jacta);
  FitnessI lafi(jacta);
  
  //Comportamiento de los modulos de los gaps se forma una "diagonal"
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

  int *s;
  s = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       if ((i%2)==0)
            s[i] = 1;
       else
            s[i] = -1;
  }
  
  cap.set_params_muestreo(samplesize, nodes, stepsize, wol, modules);
  cap.set_partition(predpar);
  bas.open_ofstream(fs, "pars/"+bas.inttostring(seed)+"_pars_set1.txt");
  cap.print_params_mues(fs);
  fs.close();
  
  int nu_unos;
  GraphI old_nw(nodes, jacta, true);
  GraphI new_nw(jacta);
  
  nu_unos = cap.make_modular_nw_unos_mues(old_nw, s, nodes);
  cap.complete_conectiviness(old_nw, s, nodes, nu_unos);
  
  

  int *mgood;
  int *mbad;
  int **oc;    
  int **ot;
  
  mgood = new int [nodes];
  mbad = new int [nodes];
  bas.create_array(oc, nodes, modules);
  bas.create_array(ot, nodes, modules);
  
  double **matmod;
  bas.create_array(matmod, nodes, nodes);
  double **prom_adjacency;
  bas.create_array(prom_adjacency, nodes, nodes);
  bas.fillmat0(prom_adjacency, nodes, nodes);
  int **temp_adjacency;
  bas.create_array(temp_adjacency, nodes, nodes);
  double modularidad;
  int edges;
  
  ofstream data;
  bas.open_ofstream(data, "hist/" + bas.inttostring(seed) + "_data_set1.txt");
  data << "Mod" << "\t" << "Edges" << endl;
  
  bas.open_ofstream(omuestra, "muestrario/" + bas.inttostring(seed) + "_muestra_set1.txt");
  bas.open_ofstream(oindex, "muestrario/" + bas.inttostring(seed) + "_inventario_set1.txt");
   for (int i = 0; i < samplesize; i++){
       if(i%(samplesize/10) == 0){
        cout << "nw " <<  i << endl;
       }
       old_nw.interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);

       do{
            cap.muestreo_rn_v2_onenw_gaps(old_nw, s, nodes, nu_unos, new_nw, mgood, mbad, oc, ot, intmods, nugaps);
            edges = new_nw.number_of_edges();
            old_nw.clear();
            old_nw.copy(new_nw);
            new_nw.clear();
       }
       while(edges != conectivity*nodes*nodes);
    
       old_nw.build_moma_d(matmod);
       modularidad = old_nw.eval_mod(matmod, predpar);
       old_nw.copy_adjacency_matrix(temp_adjacency);
       for(int i = 0; i < nodes; i++)
           for(int j = 0; j < nodes; j++)
               prom_adjacency[i][j] += temp_adjacency[i][j];
       edges = old_nw.number_of_edges();

       data << modularidad << "\t" << edges << endl;
       
       old_nw.export_nw_fixed_width_bignw(omuestra, oindex, i);
   }
   
  ofstream adj;
  bas.open_ofstream(adj, "hist/" + bas.inttostring(seed) + "_adj_set1.txt");
  for(int i = 0; i < nodes; i++){
        for(int j = 0; j < nodes; j++){
            prom_adjacency[i][j] /= samplesize;
            adj << prom_adjacency[i][j] << "\t";
        }
        adj << endl;
   }

  omuestra.close();
  oindex.close();

  data.close();
  
  for(int i = 0; i < nodes; i++){
      for(int j = 0; j < nugaps; j++)
          delete intmods[i][j];
      delete intmods[i];
  }
  delete[] intmods;
  
  delete[] mgood;
  delete[] mbad;
  for(int n = 0; n < nodes; n++){
      delete oc[n];
      delete ot[n];
  }
  delete[] oc;
  delete[] ot;

  
  jacta.close_rng();
  cout << 1 << endl;
  return 0;
}
 
