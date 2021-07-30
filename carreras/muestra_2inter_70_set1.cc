#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphi.h"
#include "muestreoi.h"
#include "fitnessi.h"


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using namespace std;

int main()
{
  
  int seed = 24125;
  int nodes = 24;
  int modules = 4;
  double conectividad = 3*nodes;
  //initial mues se refiere al número de redes en cada archivo
  int initial_mues = 50000;
  int stepsize = 1000;
  int nt_per_inter = 300;
  int inter = 2;
  int total = inter*nt_per_inter;

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI **red;
  red = new GraphI*[inter];
  for (int i = 0; i < inter; i++)
      red[i] = new GraphI[nt_per_inter];
  MuestreoI chw(jacta);
  FitnessI lafi(jacta);

  ofstream fs;
  
  chw.set_params_muestreo(total, nodes, stepsize, conectividad, modules); 
  bas.open_ofstream(fs, "pars/"+bas.inttostring(seed)+"_pars_set1.txt");
  chw.print_params_mues(fs);
  fs.close();
//   Abre los archivos de muestra
  int nu_sfiles = 2;
  string directorio[nu_sfiles];
  directorio[0] = "../premuestra";
  directorio[1] = "../premuestra/muestrario";
  

//   Exporta los parametros de muestra
  int **claves;
  bas.create_array(claves, inter, nt_per_inter);

  
  ifstream cl_anc;
  int *edges;
  edges = new int[initial_mues];
  double* modpred;
  modpred = new double[initial_mues];
  double* sort_pred;
  sort_pred = new double[initial_mues];
  int* ind_pred;
  ind_pred = new int[initial_mues];
  int count;
  

  bas.open_ifstream(cl_anc, directorio[0] + "/hist/" + bas.inttostring(seed) + "_data_set1.txt");
  count = 0;
  string header;
  cl_anc >> header >> header;
  while(cl_anc >> modpred[count] && count < initial_mues){
      cl_anc >> edges[count];
      count++;
  }
  cl_anc.close();
  bas.sort_with_index(modpred, sort_pred, initial_mues, ind_pred);
  
//   for(int i = 0; i < initial_mues; i++){
//       cout << edges[i] << endl;
//   }

        
  count = 0;
  
  for(int i = 0; i < nt_per_inter; i++){
      claves[0][i] = ind_pred[count + (int)(0.15*initial_mues)];
      count++;
  }
  
  count = 0;
  for(int i = 0; i < nt_per_inter; i++){
      claves[1][i] = ind_pred[(int)(0.85*initial_mues) - count];
      count++;
  }
    
  ofstream mod_int;
  bas.open_ofstream(mod_int, "pars/"+bas.inttostring(seed)+"_mod_int_set1.txt");
  mod_int << "1 intervalo: " << modpred[claves[0][0]] << "\t" << modpred[claves[0][nt_per_inter-1]] << endl;
  mod_int << "2 intervalo: " << modpred[claves[1][nt_per_inter-1]] << "\t" << modpred[claves[1][0]] << endl;
  mod_int.close();
      
    
  ofstream omuestra;
  ofstream oindex;
  
  bas.open_ofstream(omuestra, "muestrario/" + bas.inttostring(seed) + "_muestra_set1.txt");
  bas.open_ofstream(oindex, "muestrario/" + bas.inttostring(seed) + "_inventario_set1.txt");
  
  ifstream muestra;
  ifstream inventario;
  

  //Obtiene los genotipos extremos
  for(int j = 0; j < inter; j++){
      for (int i = 0; i < nt_per_inter; i++){
        bas.open_ifstream(muestra, directorio[1] + "/" + bas.inttostring(seed) + "_muestra_set1.txt");
        bas.open_ifstream(inventario, directorio[1] + "/" + bas.inttostring(seed) + "_inventario_set1.txt");

        red[j][i] = jacta;
        red[j][i].get_dir_nw_from_file_bignw(nodes, muestra, inventario, claves[j][i], initial_mues);
        red[j][i].export_nw_fixed_width_bignw(omuestra, oindex, i +(j*nt_per_inter));

        muestra.close();
        inventario.close();
    }
  }
  omuestra.close();
  oindex.close();
  
  set<set<int> > predpar;
  set<int>* i_par;
  i_par = new set<int>[modules];
  
  for(int j = 0; j < modules; j++){
    for (int i = j*(nodes/modules); i < (j+1)*(nodes/modules); i++)
        i_par[j].insert(i);
    predpar.insert(i_par[j]);
    i_par[j].clear();
  }
  
  double **matmod;
  bas.create_array(matmod, nodes, nodes);
  double *modularidad;
  modularidad = new double [total]; 
  int *edges_f;
  edges_f = new int[total];
  
  for(int j = 0; j < inter; j++){
      for (int i = 0; i < nt_per_inter; i++){
          red[j][i].build_moma_d(matmod);
          modularidad[(j*nt_per_inter)+i] = red[j][i].eval_mod(matmod, predpar);
          edges_f[(j*nt_per_inter)+i] = red[j][i].number_of_edges();
      }
  }
  predpar.clear();
  
  ofstream data;
  bas.open_ofstream(data, "hist/" + bas.inttostring(seed) + "_data_set1.txt");
  data << "Mod" << "\t" << "Edges" << endl;
    
  for(int j = 0; j < total; j++)
      data << modularidad[j] << "\t" << edges_f[j] << endl;
    
  data.close();
  
  jacta.close_rng();
  cout << 1 << endl;  
}
