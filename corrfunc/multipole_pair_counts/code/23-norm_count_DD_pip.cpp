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
  cout<<"Reading Input"<<endl;

  string data_file, bw_file, output_file;
  cin>>data_file;
  cin>>bw_file;
  cin>>output_file;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////

  cout<<"Loading data cat"<<endl;

  lrg gal;
  sfname<<data_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal.set_fname (fname);
  gal.fill_cat_pip();
  
  cout<<"Loading bw weights"<<endl;
 
  int cntg = (int)(gal.id.size());

  int nx = 31;
  int ny = 60;


  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);
  cout<<"bw_file: "<<bw_file<<endl;

  load_bw(bw_file, ny, bw_ids, bw_weights);



  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  double DD_norm = 0.;
  int i;
  cout<<"Calculating DD norm"<<endl;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    if (gal.clustering[i]==1) {
      for (int j=i+1;j<cntg;j++){
        if (gal.clustering[j]==1) {
          int sum_fbr = 0;
          int sum_cov = nx*ny;
          for (int is=0;is<ny;is++) {
#pragma omp atomic
            sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is) & bw_weights(j, is)));
          }
          double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
#pragma omp atomic
          DD_norm += gal.syst[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight;
        }
      }
    }
    if ((i+1)%10000==0) {
      cout<<"Completed "<<i<<" galaxies";
      print_time();
    }
  }
  cout<<"Finished DD ";
  print_time();
  
  cout<<"DD_norm: "<<DD_norm<<endl;    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"DD_norm"<<endl<<setprecision(6)<<fixed<<DD_norm;
  pfw.close();
  
  return (0);
}
