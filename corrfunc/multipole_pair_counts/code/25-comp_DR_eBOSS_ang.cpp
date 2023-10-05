#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <bitset>

#include "classes/functions.h"
#include "classes/catalogue.h"









int main (void){

  
  ostringstream sfname;
  string fname;
  ifstream pfr;
  ofstream pfw;

 
  ////////////////////////////////////////////////////////////////
  //////////          SELECT FIELD: NGC or SGC          //////////
  ////////////////////////////////////////////////////////////////
  string data_input, rand_input, bw_input, par_output, fib_output;
  cin>>data_input;
  cin>>rand_input;
  cin>>bw_input;
  cin>>par_output;
  cin>>fib_output;

  cout<<"Parameter Loading Complete"<<endl;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  lrg gal;
  sfname<<data_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal.set_fname (fname);
  gal.fill_cat_pip();
    
  int cntg = (int)(gal.id.size());

  int nx = 31;
  int ny = 60;

  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);

  load_bw(bw_input, ny, bw_ids, bw_weights);

  cout<<"Catalogue Loading Complete"<<endl;

  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rnd;
  sfname<<rand_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  rnd.set_fname (fname);
  rnd.fill_cat();

  int cntr = (int)(rnd.ra.size());

  cout<<"Random Catalogue Loaded"<<endl;


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double tmin = 0.01;
  double tmax = 3.;
  double st = 0.075;
  tmin = tmin * pow(10., -st/2.);
  int nBint = (int)(1./st*log10(tmax/tmin))+1;
  tmax = tmin * pow(10., nBint*st);
  cout<<"Binning Complete"<<endl;
  

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  vector<double> DR_par(nBint);
  vector<double> DR_fib(nBint);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++) {
    double iip_weight;
    if (gal.fiber[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
#pragma omp atomic
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i,is)));
      }
      iip_weight = (double)(sum_cov)/(double)(sum_fbr);
    }
    for (int j=0;j<cntr;j++) {
      double pair_sep = ang_sep(gal.ra[i], gal.dec[i], rnd.ra[j], rnd.dec[j]);
      if (pair_sep>tmin && pair_sep<tmax) {
        int bint = (int)(log10(pair_sep/tmin)/st);
#pragma omp atomic
        DR_par[bint] += 1.;
        if (gal.fiber[i]==1) {
#pragma omp atomic
          DR_fib[bint] += iip_weight;
        }
      }
    }
    if ((i%10000)==0){
      cout<<"I = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntg;
      print_time();
    }
  }
    
  cout<<"Pair Counting Complete"<<endl;

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<par_output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  for (int i=0;i<nBint;i++) {
    pfw<<tmin*pow(10., ((double)(i) + 0.5)*st)<<"\t"<<setprecision(6)<<fixed<<DR_par[i]<<endl;
  }
  pfw.close();

  sfname<<fib_output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  for (int i=0;i<nBint;i++) {
    pfw<<tmin*pow(10., ((double)(i) + 0.5)*st)<<"\t"<<setprecision(6)<<fixed<<DR_fib[i]<<endl;
  }
  pfw.close();

  
  cout<<"Saving Complete"<<endl;

  return (0);
}
