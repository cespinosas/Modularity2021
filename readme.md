# On the effects of the modularity of gene regulatory networks on phenotypic variability and its association with robustness
Ulises Hernández, Luis Posadas-Vidales, and Carlos Espinosa-Soto*

Instituto de Fı́sica, Universidad Autónoma de San Luis Potosı́, Manuel Nava 6, Zona 
Universitaria, San Luis Potosı́, Mexico
May 19, 2021

## Description of contents
Here you can find all the code required to produce the results reported in the main text and supplementary
material of the article “On the effects of the modularity of gene regulatory networks on phenotypic variability
and its association with robustness”, authored by us. All scripts are written in C++ and R. The programs
require the GSL numeric library for C++ (https://www.gnu.org/software/gsl/). Therefore,
compilation using g++ requires including the flags -lgsl -lgslcblas.
Our codes saves the data that we generate in a specific directory structure that we provide here.

## libs/
This directory contains libraries with functions that our scripts require. Compilation of most of our code
requires linking with some of the libraries in this directory.

## premuestra/
This directory contains the programs required to build large samples of networks. The muestrario/
directory is used to store the networks; the hist/ directory is used to store data related to the sample; and
the pars/ directory is required to save the parameters used to build the sample. The premuestra_\*.cc
files provide the code required to build the networks. Hereafter, if the term “set1” is in the name of a file,
the networks on that file sustain the target GAPs. On the contrary, if the term “2m” is in the name of a file,
the networks on that file sustain the alternative target GAPs. The \*.R files are required to plot graphs of
properties of the sample.

## robustez\_miope/
This directory contains the programs required to study the association between the modularity and the
robustness of a network. Moreover, the programs in this directory allow us to study the GAPs produced by the
single-mutants of a network. The name of the subdirectories within this directory refers to the experiments
performed within each subdirectory. The access_phen/ subdirectory contains the programs required
to study the GAPs of the single-mutants of a network. The miope/ subdirectory provides the programs
required to perform the directed evolution of a network to increase their robustness or modularity. The
robustez/ subdirectory provides the programs required to study the correlation between the modularity
and robustness of a network. Within each directory the muestra_\*.cc files contain the code required to
move already sampled networks (in the premuestra/ directory) into the appropriate location.
Within the access_phen/ subdirectory, the 1mut_access_phen*.cc files are required to build
the GAPs of the single-mutants of a network and analyze their characteristics. These files use the access_phen/
directory to export their results. Within the miope/ subdirectory, the miope_mod_all_*.cc files
provides the code that is used to evolve a network under selection for greater modularity. On the other hand,
the miope_rob_all_*.cc files provides the code required to evolve a network under selection for greater
robustness. The files export their resulting data into the miope_mod/ directory or into the miope_rob/
subdirectory depending on the criteria used to select the networks during their evolution: modularity or
robustness respectively. Within the robustez/ subdirectory, the rob_dist_1mut_multic_\*.cc
files contain the code required to assess the correlation between the modularity of a network and its modularity.
The files export their results into the rob/ directory. Within each subdirectory the \*.R files contain the code
required to graph the results produced by each experiment.

## modab
This directory contains tools that allow us to analyze the effect of modularity in the phenotypic variability of
networks under the effect of the accumulation of multiple mutations. The muestra_2inter_70_\*.cc
files contain code required to build networks sets with high and with low modularity using the sample in
the premuestra directory. These files export the selected set of networks to the modab/premuestra/
subdirectory. The file mult_carp_\*.cc contains the code to prepare the aforementioned set of networks
for the following analyses. The distbtnw_multrep_transph_\*.cc files provide the code required
to analyze the effect of modularity in the phenotypic variability of a network under the accumulation of
multiple mutations. The \*.R files are used to graph the results produced by the latter programs

## carreras
This directory contains the programs required to study simulations of the evolution of populations of GRNs. Within 
this directory the \*.cc files contains the code required to select networks within a modularity interval from 
already sampled networks (in the premuestra/ directory) and move them into the appropriate location. The
carreras\_m1nuic/ subdirectory contains three folders. The onepeak\_1chan/ folder contains code required to 
perform evolutionary simulations of populations of GRN where the optimal GAP has one change with respect to 
the ancestral GAP. The onepeak_2nomodchanges contains code required to perform evolutionary simulations of populations 
of GRN where the optimal GAP has two changes with respect to the ancestral GAP. Each of these changes is located in a
different coordinated group. The onepeak_2modchanges contains code required to perform evolutionary simulations of
populations of GRN where the optimal GAP has two changes with respect to the ancestral GAP. Both changes are located
in the same coordinated group. Within each of these three directories there are four program files. The previo.cc file
provides the code to generate the folders required during evolutionary simulations. The digievol.cc file provides
the code to perform the evolutionary simulation of a population of GRNs. This program requires as an argument an
integer to identify the GRN that will be used as founder for the population. Such founder is stored in the 
carreras/muestrario/ subdirectory. The previo.cc file and the time_which.R file contains the code required to analyze
the number of generations that the populations take to found the optimal GAP.

