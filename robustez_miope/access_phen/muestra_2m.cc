#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
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
using std::string;
using namespace std;

int main()
{
  
  int seed = 24125;
  int nodes = 24;
  int modules = 4;
  double conectividad = 3*nodes;
  //initial mues se refiere al número de redes en cada archivo
  int initial_mues = 5000;
  int stepsize = 20*nodes*nodes;
  int sample_nt = 5000;

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI *red;
  red = new GraphI[sample_nt];
  Cap_evolI chw(jacta);
  FitnessI lafi(jacta);

  ofstream fs;
  
  chw.set_params_muestreo(sample_nt, nodes, stepsize, conectividad, modules); 
  bas.open_ofstream(fs, "pars/"+bas.inttostring(seed)+"_pars_2m.txt");
  chw.print_params_mues(fs);
  fs.close();
//   Abre los archivos de muestra
  int nu_sfiles = 2;
  string directorio[nu_sfiles];
  directorio[0] = "../../premuestra";
  directorio[1] = "../../premuestra/muestrario";
    
  ofstream omuestra;
  ofstream oindex;
  
  bas.open_ofstream(omuestra, "muestrario/" + bas.inttostring(seed) + "_muestra_2m.txt");
  bas.open_ofstream(oindex, "muestrario/" + bas.inttostring(seed) + "_inventario_2m.txt");
  
  ifstream muestra;
  ifstream inventario;
  
  int warn_edges;
  int count = 0;
  //Obtiene los genotipos extremos
  for (int i = 0; i < sample_nt; i++){
      bas.open_ifstream(muestra, directorio[1] + "/" + bas.inttostring(seed) + "_muestra_2m.txt");
      bas.open_ifstream(inventario, directorio[1] + "/" + bas.inttostring(seed) + "_inventario_2m.txt");

      red[i] = jacta;
      red[i].get_dir_nw_from_file_bignw(nodes, muestra, inventario, count, initial_mues);
      warn_edges = red[i].number_of_edges();
      while(warn_edges != conectividad){
        count++;
        red[i].clear();
        red[i].get_dir_nw_from_file_bignw(nodes, muestra, inventario, count, initial_mues);
        warn_edges = red[i].number_of_edges();
//         cout << warn_edges << endl;
        if(warn_edges < conectividad - 4 || warn_edges > conectividad + 4){
            cout << "[Error]: Initial sample connectivity does not contain wanted connectivity.\n";
            exit(1);
        }
      }
      red[i].export_nw_fixed_width_bignw(omuestra, oindex, i);
      count++;
      
      muestra.close();
      inventario.close();
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
  modularidad = new double [sample_nt]; 
  int *edges_f;
  edges_f = new int[sample_nt];
  
  for (int i = 0; i < sample_nt; i++){
      red[i].build_moma_d(matmod);
      modularidad[i] = red[i].eval_mod(matmod, predpar);
      edges_f[i] = red[i].number_of_edges();
  }
  predpar.clear();
  
  ofstream data;
  bas.open_ofstream(data, "hist/" + bas.inttostring(seed) + "_data_2m.txt");
  data << "Mod" << "\t" << "Edges" << endl;
    
  for(int j = 0; j < sample_nt; j++)
      data << modularidad[j] << "\t" << edges_f[j] << endl;
    
  data.close();
  
  jacta.close_rng();
  cout << 1 << endl;  
}
