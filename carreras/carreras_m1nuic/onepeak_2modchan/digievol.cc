#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphi.h"
#include "evolvei.h"
#include "muestreoi.h"
#include "fitnessi.h"


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2) {
        cout << "[Error]; Wrong number of arguments in program. It should be ancester number \n";
        exit(1);
        
    }
  
  double multexp = 1;
  int cl_ances = atoi(argv[1]);
  int seed = cl_ances;
  int seed_ances = 24125;
  int nodes = 24;
  double wol = double(3)/double(nodes);
  int modules = 4;
  int genpermod = nodes/modules;
  
  int popsize = 1000;
  double numutgenes = 0.8;
  double selcoefpg = 0.4;
  double max_generation = 100000000;
  int chan_genes = 2;
  int nugoales = 1;
  int at_step = 200;
  int gen_mues = 1000;
  
  int inter_mod = 2;
  double initial_mb = 0.0968;
  double initial_ma = 0.1930;
  double mues_sd = 0.0464251524179295;
  
  double murate = 1.0 - pow(numutgenes, (1.0/double(nodes)));


  Alea jacta(seed);
  Basics bas(jacta);
  EvolveI chw(jacta);
  GraphI red(jacta);
  GraphC redprom(jacta);
  GraphC redprom_1(jacta);
  GraphI redmore(jacta);
  GraphI redmore_1(jacta);
  FitnessI fifu(jacta);
  MuestreoI mues(jacta);
    
  //Abre los archivos de muestra
  int nu_sfiles = 7;
  string directorio[nu_sfiles];
  directorio[0] = "../../muestrario";
  directorio[1] = "primer_opt";
  directorio[2] = "25_opt";
  directorio[3] = "wtrack";
  directorio[4] = "dead";
  directorio[5] = "trunc";
  directorio[6] = "gen_dist";
   
  //produce the two genes that have to change and the gaps
  //the next code functions if only one ic is evaluated. We don't care for the rest of the gaps. 
  //*note: an alternative is asking the networks to produce the m gaps (m=modules) and aditionally find randomly another gap. This gap would be the one we select for.
  
  int nugaps = modules;
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
  
  int *s;
  s = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       if ((i%2)==0)
            s[i] = 1;
       else
            s[i] = -1;
  }
  
  int **all_ats;
  bas.create_array(all_ats, nugaps, nodes);
  mues.start_gaps_ics(all_ats, nugaps, modules, nodes, modgaps);
  
  int *mgood, *mbad;
  int **oc, **ot;
  
  mgood = new int[nodes];
  mbad = new int[nodes];
  oc = new int*[nodes];
  ot = new int*[nodes];
  for(int i = 0; i < nodes; i++){
      oc[i] = new int[modules];
      ot[i] = new int[modules];
  }
  
  //lee la red ancestro
  ifstream muestra;
  bas.open_ifstream(muestra, directorio[0] + "/" + bas.inttostring(cl_ances) + "/nw_set1.txt");
  red.get_dir_nw_from_file(nodes, muestra);
  muestra.close();
  
  redprom_1.copy(red);
  redmore_1.copy(red);
  
  int *ic;
  ic = new int[nodes];
  for (int i = 0; i < nodes; i++) {
       ic[i] = all_ats[(modules-1)][i];
  }

  int** genes_co;
  bas.create_array(genes_co, nugoales, chan_genes);
  genes_co[0][0] = 9;
  genes_co[0][1] = 10;
  
  int **goals;
  bas.create_array(goals, nugoales, nodes);
  for (int i = 0; i < nodes; i++) {
          goals[0][i] = ic[i];
    }
  goals[0][genes_co[0][0]] *= -1;  
  goals[0][genes_co[0][1]] *= -1;
  
  //fija las modularidades
  set<set<int> > predpar;
  set<int>* i_par;
  i_par = new set<int>[modules];
  
  for(int j = 0; j < modules; j++){
    for (int i = j*(nodes/modules); i < (j+1)*(nodes/modules); i++)
        i_par[j].insert(i);
    predpar.insert(i_par[j]);
    i_par[j].clear();
  }
  delete [] i_par;
  
  double *to_mod, *from_mod;
  to_mod = new double[inter_mod];
  from_mod = new double[inter_mod];
  to_mod[0] = initial_mb + mues_sd;
  from_mod[0] = initial_mb - mues_sd;
  to_mod[1] = initial_ma + mues_sd;
  from_mod[1] = initial_ma - mues_sd;
  
  double **matmod;
  double mod, modmin, modmax;
  bas.create_array(matmod, nodes, nodes);
  red.build_moma_d(matmod);
  mod = red.eval_mod(matmod, predpar);
  cout << mod << endl;
  if(mod <= to_mod[0] && mod >= from_mod[0]){
      modmin = from_mod[0];
      modmax = to_mod[0];
  }  
  else if(mod <= to_mod[1] && mod >= from_mod[1]){
      modmin = from_mod[1];
      modmax = to_mod[1];
  }
  else{
      cout << "[Error]: digievol ancester network has a modularity not specified in intervals\n";
      exit(1);
  }
  
  ofstream fs;
  chw.set_params(murate, wol, popsize, nodes);      
  bas.open_ofstream(fs, "pars/"+bas.inttostring(seed_ances)+"_pars.txt");
  chw.print_params(fs);
  fs.close();
  ofstream opf1t;
  ofstream eopn1;
  ofstream fsw, fsnoter, fsd, fst, fgd;
  
  bas.open_ofstream(eopn1, directorio[1] + "/redes/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + ".txt");
  
  bas.open_ofstream(opf1t, directorio[1] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances) + "_time.txt");

  bas.open_ofstream(fsw, directorio[3] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_w.txt");

  bas.open_ofstream(fsnoter, directorio[4] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_noter.txt");
  
  bas.open_ofstream(fsd, directorio[4] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_dead.txt");
  fst << "muertos_emb"  << "\t" << "muertos_ad" << endl;

  bas.open_ofstream(fst, directorio[5] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_trunc.txt");
  fst << "individuo" << "\t" << "generacion" << endl;
  
  bas.open_ofstream(fgd, directorio[6] + "/" + bas.inttostring(seed_ances) + "_" + bas.inttostring(cl_ances)  + "_dist.txt");
  fgd << "Distancia_ances_mean" << "\t" << "Distancia_ances_more" << "\t" << "Porcentaje_more" << "\t" << "Velocidad_mean" << "\t" << "Velocidad_more" << "\t" << "Desparramamiento" << endl;
  
  set<int> losopt;
  set<int>::iterator it;  
        
  chw = jacta;
  chw.set_params(murate, wol, popsize, nodes);      
  chw.start_pop(red);
    
    
  double w, wpre;
  double d;
  double llego = -1;
  int count_ge = 0;
  double maxw = 0;
  int vivio = 0;
  int dead_em = 0;
  int dead_ad = 0;
  double porc_more, desparrame;

  do {
      
      for (int j = 0; j < popsize; j++){
            //condición de mod
            chw.population[j].build_moma_d(matmod);
            mod = chw.population[j].eval_mod(matmod, predpar);
            w = 0;
            wpre = 0;
            chw.population[j].interactions_to_attractor(s, mgood, mbad, oc, ot, nodes, modules, predpar);
            if(fifu.cap_to_s_gaps(nodes, mgood, mbad, oc, ot, modules, intmods, (nugaps-1))){
                if(mod < modmax && mod > modmin){
                    chw.population[j].set_as_state(ic);
                    if(chw.population[j].find_an_attractor_trunc(at_step)){
                        d = 1.0 - chw.population[j].distance_pub(goals[0]);
                        d *= (multexp*nodes);
                        wpre = pow((1-selcoefpg), d);
                    }
                    else{
                        fst  << j << "\t" << count_ge << endl;
                    }
                    if(wpre == 1)
                        w = 1;
                    else{
                        d = 1.0 - chw.population[j].distance_pub(ic);
                        d *= (multexp*nodes);
                        d += chan_genes;
                        w = pow((1-selcoefpg), d);
                    }
                    if (w == 1){
                        if(llego == -1)
                            llego += 1;
                    }
                    if (llego == 0){
                        opf1t << count_ge << endl;
                        chw.population[j].export_nw(eopn1);                
                        llego += 1;
                    }
                    chw.population[j].clear_attractor();
                }
            }
            dead_ad++;
            chw.assign_w(j, w);
      }   
        
      chw.calc_meanw();
      maxw = chw.return_maxw();
      fsw << chw.return_meanw() << "\t" << maxw << "\t" << chw.return_sdw() << endl;

      //Muestreos de redes en la población
      if(count_ge % gen_mues == 0){
          redprom.make_nw(nodes, true);
          redmore.make_nw(nodes, true);
          chw.mean_nw(redprom);
          porc_more = chw.more_repeated_nw(redmore);
          desparrame = chw.dist_bw_nws();
          
          fgd << redprom.distance_between_nw(red) << "\t" << redmore.distance_between_nw(red) << "\t" << porc_more << "\t" << redprom_1.distance_between_nw(redprom) << "\t" << redmore_1.distance_between_nw(redmore) << "\t" << desparrame << endl;
          
          redprom_1.clear();
          redmore_1.clear();
          redprom_1.copy(redprom);
          redmore_1.copy(redmore);
          redprom.clear();
          redmore.clear();
      }
    
      if (count_ge < max_generation && llego != 1 && maxw != 0)
          dead_em += chw.one_generation_mod_dead(modmax, modmin, predpar);
      count_ge++;
      vivio++;
  } while(count_ge <= max_generation && llego != 1 && maxw != 0);
    
  if (maxw == 0){
      opf1t << "dead" << endl;
      fsnoter << vivio;
  }
  else if(count_ge > max_generation){
      opf1t << "not_enough_time" << endl;
      fsnoter << "not_enough_time";
  }
  else 
      fsnoter << "NA" << endl;

  fsd << dead_em << "\t" << dead_ad << endl;
  
  fsw << endl;
  
  fsw.close();
  fsnoter.close();
  fsd.close();
  fst.close();
  fgd.close();

  red.clear();
  jacta.close_rng();
  cout << 1 << endl;
  return 0;

}
