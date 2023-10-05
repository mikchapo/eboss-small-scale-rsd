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
  string data_file, rand_file, bw_file, output_file;
  cin>>data_file;
  cin>>rand_file;
  cin>>bw_file;
  cin>>output_file;
  

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  lrg gal;
  sfname<<data_file;
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
  cout<<"bw_file: "<<bw_file<<endl;

  load_bw(bw_file, ny, bw_ids, bw_weights);

  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rnd;
  sfname<<rand_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  rnd.set_fname(fname);
  rnd.fill_cat();

  int cntr = (int)(rnd.ra.size());


  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
//  double DR_norm = 0.;
//  int i,j;   
//#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
//  for (i=0;i<cntg;i++){
//    if (gal.clustering[i]==1) {
//      for (j=0;j<cntr;j++){
//        int sum_fbr = 0;
//        int sum_cov = nx*ny;
//        for (int is=0;is<ny;is++) {
//#pragma omp atomic
//          sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
//        }
//        double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
//#pragma omp atomic
//        DR_norm += gal.syst[i]*gal.noz[i]*gal.fkp[i]*rnd.syst[j]*rnd.cp[j]*rnd.noz[j]*rnd.fkp[j]*pip_weight;
//      }
//    }
//  }
//  cout<<"Finished DR ";
//  print_time();

  double D_eff = 0.;
  double R_eff = 0.;
  int i;   
  for (i=0;i<cntg;i++){
    if (gal.clustering[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
      }
      double iip_weight = (double)(sum_cov)/(double)(sum_fbr);
      D_eff += gal.syst[i]*gal.noz[i]*gal.fkp[i]*iip_weight;
    }
  }
  for (i=0;i<cntr;i++){
    R_eff += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
  }
  cout<<"Finished DR ";
  print_time();

  double DR_norm = D_eff*R_eff;
  cout<<"D_eff: "<<D_eff<<" R_eff: "<<R_eff<<" DR_norm: "<<DR_norm<<endl;
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"D_eff\tR_eff\tDR_norm"<<endl<<setprecision(6)<<fixed<<D_eff<<"\t"<<R_eff<<"\t"<<DR_norm;
  pfw.close();
  
  return (0);
}
