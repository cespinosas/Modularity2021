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
  int muestras = 200;

  Alea jacta(seed);
  Basics bas(jacta);
  GraphI red1(jacta);
  GraphI red2(jacta);
  Cap_evolI chw(jacta);

  string directorio = "muestrario";
  

  for (int i = 0; i < muestras; i++)
    system(("mkdir -p " + directorio + "/" + bas.inttostring(i)).c_str());

    
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
//   Exporta los parametros de muestra
  ifstream muestra;
  ifstream inventario;
  ofstream omuestra;

  //Obtiene los genotipos extremos
  for (int i = 0; i < muestras; i++){
    bas.open_ifstream(muestra, directorio + "/" + bas.inttostring(seed) + "_muestra_2m.txt");
    bas.open_ifstream(inventario, directorio + "/" + bas.inttostring(seed) + "_inventario_2m.txt");

    red1.get_dir_nw_from_file_bignw(nodes, muestra, inventario, i, muestras);
    
    muestra.close();
    inventario.close();
    
    red1.build_moma_d(matmod);
    cout << red1.eval_mod(matmod, predpar) << endl;
          
    bas.open_ofstream(omuestra, directorio + "/" + bas.inttostring(i) + "/" + "nw_2m.txt");
    red1.export_nw(omuestra);

    omuestra.close();
    red1.clear();
  }

  jacta.close_rng();
  cout << 1 << endl;  
}
