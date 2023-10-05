#ifndef CATALOGUE_H
#define CATALOGUE_H






#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <cmath>
#include <iomanip>

#include <algorithm>






using namespace std;




class catalogue{
 protected:
  string fname;
 public:
  void set_fname(const string fnamein){fname = fnamein;}
};





class lrg: public catalogue{
 public:
  vector<long int> id;
  vector<int> fiber, clustering;
  vector<double> ra,dec,zs,syst,cp,fkp,noz,x,y,z,jk;
  void fill_cat_basic();
  void fill_cat_ang();
  void fill_cat_ang_weighted();
  void fill_cat();
  void fill_cat_weighted();
  void fill_cat_jk();
  void fill_cat_pip();
};



class randcat: public catalogue{
 public:
  vector<double> ra,dec,zs,syst,cp,fkp,noz,x,y,z,jk;
  void fill_cat_basic();
  void fill_cat_ang();
  void fill_cat();
  void fill_cat_jk();
};


class weights{
 protected:
  string fname;
 public:
  void set_fname(const string fnamein){fname = fnamein;}
};


class ang_weights: public weights{
 public:
  vector<double> theta,ang;
  void load_weights_basic();
  void load_weights();
};


#endif
