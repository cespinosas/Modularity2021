#ifndef BASICS_H
#define BASICS_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include <bitset>
#include "alea.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <sys/stat.h>
#include <dirent.h>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

class Basics
{
  public:
  Basics();
  Basics(Alea& jacta);
  void start_rng(Alea& jacta);

  //create array
  void create_array(int** &arr, int rows, int cols);
  void create_array(bool** &arr, int rows, int cols);
  void create_array(double** &arr, int rows, int cols);
  void create_array(char** &arr, int rows, int cols);
  void create_array(string** &arr, int rows, int cols);
  void create_array(int*** &arr, int slices, int rows, int cols);
  void create_array(bool*** &arr, int slices, int rows, int cols);
  void create_array(double*** &arr, int slices, int rows, int cols);
  void create_array(char*** &arr, int slices, int rows, int cols);
  void create_array(string*** &arr, int slices, int rows, int cols);
  void create_array(int**** &arr, int boxes, int slices, int rows, int cols);
  void create_array(double**** &arr, int boxes, int slices, int rows, int cols);
  void create_array(bool**** &arr, int boxes, int slices, int rows, int cols);
  void create_array(int***** &arr, int rooms, int boxes, int slices, int rows, int cols);

  //open-close files
  void open_ifstream(ifstream& fe, string nomb);
  void open_ofstream(ofstream& fs, string nomb);
  void open_ofstream_to_append(ofstream& fs, string nomb);

  //vectors
  void fillv0(int vec[], int s);
  void fillv0(double vec[], int s);
  void fillv0(bool vec[], int s);
  void fillv0(string vec[], int s);
  void fillv1(int vec[], int s);
  void fillv1(double vec[], int s);
  void fillv1(bool vec[], int s);
  void fillvm1(int vec[], int s);
  void fillvm1(double vec[], int s);

  bool eqvec(int vec1[], int s1, int vec2[], int s2);
  bool eqvec(double vec1[], int s1, double vec2[], int s2);
  bool eqvec(bool vec1[], int s1, bool vec2[], int s2);

  int find_max(int* vec, int tam);
  double find_max(double* vec, int tam);
  int find_min(int* vec, int tam);
  double find_min(double* vec, int tam);
  int find_max_index(int* vec, int tam);
  int find_max_index(double* vec, int tam);
  int find_min_index(int* vec, int tam);
  int find_min_index(double* vec, int tam);
  void sort(int* vec, int* nvec, int tam);
  void sort_with_index(int* vec, int* nvec, int tam, int* index);
  void sort(double* vec, double* nvec, int tam);
  void sort_with_index(double* vec, double* nvec, int tam, int* index);

  int sumatoria(int vec[], int s);
  double sumatoria(double vec[], int s);
  double get_mean(int* vec, int tam);
  double get_mean(double* vec, int tam);
  double get_sample_variance(int* vec, int tam);
  double get_sample_variance(double* vec, int tam);
  double get_pop_variance(int* vec, int tam);
  double get_pop_variance(double* vec, int tam);
  double get_sample_stddev(int* vec, int tam);
  double get_sample_stddev(double* vec, int tam);
  double get_pop_stddev(int* vec, int tam);
  double get_pop_stddev(double* vec, int tam);
  double get_sample_stderr(int* vec, int tam);
  double get_sample_stderr(double* vec, int tam);


  //matrix
  void fillmat0(int** mat, int rows, int cols);
  void fillmat0(double** mat, int rows, int cols);
  void fillmat0(bool** mat, int rows, int cols);
  void fillmatm1(int** mat, int rows, int cols);
  void fillmatm1(double** mat, int rows, int cols);
  int find_max(int** mat, int rows, int cols);
  double find_max(double** mat, int rows, int cols);

  //vector-matrix
  int vecinmat(int** mat, int rows, int cols, int vec[], int s);
  int vecinmat(double** mat, int rows, int cols, double vec[], int s);
  int vecinmat(bool** mat, int rows, int cols, bool vec[], int s);

  //strings
  string inttostring(int num);
  string bintostring(int num, int spaces);
  
  //networks
  void get_dot_fig(string arch);
  void get_dot_fig(string arch, string ext);
  void get_circo_dot_fig(string arch);
  void get_circo_dot_fig(string arch, string ext);
  void polar_to_cartesian(double radians, double radius, double& x, double& y);



  private:
  Alea est;
};

#endif



