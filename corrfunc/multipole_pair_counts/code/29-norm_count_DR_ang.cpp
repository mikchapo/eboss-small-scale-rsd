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
  double DR_norm_par = (double)(cntg)*(double)(cntr);
  double D_eff = 0.;
  double R_eff = (double)(cntr);
  int i;   
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    if (gal.fiber[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
#pragma omp atomic
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
      }
//    if (sum_fbr==0) {
//      cout<<"Object "<<gal.id[i]<<" sum_fbr=0"<<endl;
//    }
      double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
#pragma omp atomic
      D_eff += pip_weight;
    }
  }
  cout<<"Finished DR ";
  print_time();

  double DR_norm_fib = D_eff*R_eff;

  cout<<"DR_norm_par: "<<DR_norm_par<<endl<<"DR_norm_fib: "<<DR_norm_fib<<endl;
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"DR_norm_par\tDR_norm_fib"<<endl<<setprecision(6)<<fixed<<DR_norm_par<<"\t"<<DR_norm_fib;
  pfw.close();
  
  return (0);
}
