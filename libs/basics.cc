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
#include "basics.h"
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


Basics::Basics()
{
}

Basics::Basics(Alea& jacta)
{
	start_rng(jacta);
}

void Basics::start_rng(Alea& jacta)
{
  est = jacta;
}

void Basics::create_array(int** &arr, int rows, int cols) {
  int i;
  arr = new int*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new int[cols];
  return;
}

void Basics::create_array(bool** &arr, int rows, int cols) {
  int i;
  arr = new bool*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new bool[cols];
  return;
}

void Basics::create_array(double** &arr, int rows, int cols) {
  int i;
  arr = new double*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new double[cols];
  return;
}

void Basics::create_array(char** &arr, int rows, int cols) {
  int i;
  arr = new char*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new char[cols];
  return;
}

void Basics::create_array(string** &arr, int rows, int cols) {
  int i;
  arr = new string*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new string[cols];
  return;
}

void Basics::create_array(int*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new int**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new int*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new int[cols];
  }
  return;
}

void Basics::create_array(bool*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new bool**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new bool*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new bool[cols];
  }
  return;
}

void Basics::create_array(double*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new double**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new double*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new double[cols];
  }
  return;
}

void Basics::create_array(char*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new char**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new char*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new char[cols];
  }
  return;
}

void Basics::create_array(string*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new string**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new string*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new string[cols];
  }
  return;
}

void Basics::create_array(int**** &arr, int boxes, int slices, int rows, int cols) {
  int n, i, j;
  arr = new int*** [boxes];
  for (n  = 0; n < boxes; n++){
      arr[n] = new int** [slices];
    for (i = 0; i < slices; i++) {
        arr[n][i] = new int*[rows];
        for (j = 0; j < rows; j++)
        arr[n][i][j] = new int[cols];
    }
  }
  return;
}

void Basics::create_array(double**** &arr, int boxes, int slices, int rows, int cols) {
  int n, i, j;
  arr = new double*** [boxes];
  for (n  = 0; n < boxes; n++){
      arr[n] = new double** [slices];
    for (i = 0; i < slices; i++) {
        arr[n][i] = new double*[rows];
        for (j = 0; j < rows; j++)
        arr[n][i][j] = new double[cols];
    }
  }
  return;
}

void Basics::create_array(bool**** &arr, int boxes, int slices, int rows, int cols) {
  int n, i, j;
  arr = new bool*** [boxes];
  for (n  = 0; n < boxes; n++){
      arr[n] = new bool** [slices];
    for (i = 0; i < slices; i++) {
        arr[n][i] = new bool*[rows];
        for (j = 0; j < rows; j++)
        arr[n][i][j] = new bool[cols];
    }
  }
  return;
}

void Basics::create_array(int***** &arr, int rooms, int boxes, int slices, int rows, int cols) {
  int m, n, i, j;
  arr = new int**** [rooms];
  for (m = 0; m < rooms; m++){
    arr[m] = new int*** [boxes];
    for (n  = 0; n < boxes; n++){
        arr[m][n] = new int** [slices];
        for (i = 0; i < slices; i++) {
            arr[m][n][i] = new int*[rows];
            for (j = 0; j < rows; j++)
            arr[m][n][i][j] = new int[cols];
        }
    }
  }
  return;
}


void Basics::open_ifstream(ifstream& fe, string nomb)
{
	fe.open(nomb.c_str());
	if (!fe.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ifstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str());
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream_to_append(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str(), ios::app);
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

//vectors

void Basics::fillv0(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = false;
}

void Basics::fillv0(string vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = "";  
}

void Basics::fillv1(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = true;
}

void Basics::fillvm1(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

void Basics::fillvm1(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

bool Basics::eqvec(int vec1[], int s1, int vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(double vec1[], int s1, double vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(bool vec1[], int s1, bool vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

int Basics::find_max(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

double Basics::find_max(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

double Basics::find_min(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

int Basics::find_min(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

void Basics::sort(int* vec, int* nvec, int tam)
{
	int min, max;
	int i, j, quedan, k;
	int* wov;
	wov = new int[tam];
	int* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new int[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new int[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

void Basics::sort_with_index(int* vec, int* nvec, int tam, int* index)
{
	int min, max;
	int i, j, quedan, k;
	int* wov;
	wov = new int[tam];
	int* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
    
    int* inx;
	inx = new int[tam];
    int* inx2;
	for (i=0; i < tam; i++)
		inx[i] = i;
    
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
                index[i] = inx[j];
				i++;
				break;
			}
		}
		wov2 = new int[quedan-1];
        inx2 = new int[quedan-1];
		for (k = 0; k < j; k++){
			wov2[k] = wov[k];
            inx2[k] = inx[k];
        }
		for (k=j; k < (quedan-1); k++){
			wov2[k] = wov[k+1];
            inx2[k] = inx[k+1];
        }
		delete [] wov;
        delete [] inx;
		quedan--;
		wov = new int[quedan];
        inx = new int[quedan];
		for (k=0; k < quedan; k++){
			wov[k] = wov2[k];
            inx[k] = inx2[k];
        }
		delete [] wov2;
        delete [] inx2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

void Basics::sort(double* vec, double* nvec, int tam)
{
	double min, max;
	int i, j, quedan, k;
	double* wov;
	wov = new double[tam];
	double* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new double[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new double[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}


void Basics::sort_with_index(double* vec, double* nvec, int tam, int* index)
{
	double min, max;
	int i, j, quedan, k;
	double* wov;
	wov = new double[tam];
	double* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
    
    int* inx;
	inx = new int[tam];
    int* inx2;
	for (i=0; i < tam; i++)
		inx[i] = i;
    
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
                index[i] = inx[j];
				i++;
				break;
			}
		}
		wov2 = new double[quedan-1];
        inx2 = new int[quedan-1];
		for (k = 0; k < j; k++){
			wov2[k] = wov[k];
            inx2[k] = inx[k];
        }
		for (k=j; k < (quedan-1); k++){
			wov2[k] = wov[k+1];
            inx2[k] = inx[k+1];
        }
		delete [] wov;
        delete [] inx;
		quedan--;
		wov = new double[quedan];
        inx = new int[quedan];
		for (k=0; k < quedan; k++){
			wov[k] = wov2[k];
            inx[k] = inx2[k];
        }
		delete [] wov2;
        delete [] inx2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

int Basics::sumatoria(int vec[], int s)
{
	int res=0;
	int i;
	for (i=0; i<s; i++)
		res = res + vec[i];
	return res;
}

double Basics::sumatoria(double vec[], int s)
{
	double res=0;
	int i;
	for (i=0; i<s; i++)
		res = res + vec[i];
	return res;
}

double Basics::get_mean(int* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}

double Basics::get_mean(double* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}
//matrix
void Basics::fillmat0(int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(bool** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmatm1(int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

void Basics::fillmatm1(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

int Basics::find_max(int** mat, int rows, int cols)
{
  int* maxs; 
  maxs = new int[rows];
  int res;

  for (int i = 0; i < rows; i++)
      maxs[i] = find_max(mat[i], cols);

  res = find_max(maxs, rows);
    
  return res;
}

double Basics::find_max(double** mat, int rows, int cols)
{
  double* maxs;  
  maxs = new double[rows];
  double res;

  for (int i = 0; i < rows; i++)
      maxs[i] = find_max(mat[i], cols);

  res = find_max(maxs, rows);
    
  return res;
}

//vectors and matrices
int Basics::vecinmat(int** mat, int rows, int cols, int vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}


int Basics::vecinmat(double** mat, int rows, int cols, double vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols)){
				res = i;
				break;
			}
  }
  return res;
}

int Basics::vecinmat(bool** mat, int rows, int cols, bool vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s)   {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}

string Basics::inttostring(int num)
{
	char buff[50];
	int j;
	j = sprintf(buff, "%d", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to int was not possible using Basics::inttostring.\n";
		exit(1);
	}
	string res(buff);
	return res;
}

string Basics::bintostring(int num, int spaces)
{
	string spre = std::bitset< 32 >( num ).to_string();
	string s;
    
    s = spre.substr(32 - spaces);

    return s;
}
