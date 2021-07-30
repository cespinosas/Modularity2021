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
  

  int wol = (int)nodes*nodes*conectivity;
  Alea jacta(seed);
  Basics bas(jacta);
  Cap_evolI cap(jacta);
  FitnessI lafi(jacta);
  
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
  bas.open_ofstream(fs, "pars/"+bas.inttostring(seed)+"_pars_2m.txt");
  cap.print_params_mues(fs);
  fs.close();
  
  int nu_unos;
  GraphI old_nw(nodes, jacta, true);
  GraphI new_nw(nodes, jacta, true);
  
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
  bas.open_ofstream(data, "hist/" + bas.inttostring(seed) + "_data_2m.txt");
  data << "Mod" << "\t" << "Edges" << endl;
  
  bas.open_ofstream(omuestra, "muestrario/" + bas.inttostring(seed) + "_muestra_2m.txt");
  bas.open_ofstream(oindex, "muestrario/" + bas.inttostring(seed) + "_inventario_2m.txt");
   for (int i = 0; i < samplesize; i++){
       if(i%(samplesize/10) == 0){
        cout << i << endl;
       }
       old_nw.interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
//         if(!lafi.cap_to_s(nodes, mgood, mbad, oc, ot, modules)){
//             cout << "[Error]: Error in premuestra_bignw.\n";
//             exit(1);
//         }
       do{
            cap.muestreo_rn_v2_onenw(old_nw, s, nodes, nu_unos, new_nw, mgood, mbad, oc, ot);
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
  bas.open_ofstream(adj, "hist/" + bas.inttostring(seed) + "_adj_2m.txt");
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
  
  delete[] s;
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
 
