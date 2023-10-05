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

  double sepmin = 0.1;
  double sepmax = 60.;


  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  double norm = 0.;
  double z_eff = 0.;
  int i;
  cout<<"Calculating DD norm"<<endl;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    if (gal.clustering[i]==1) {
      for (int j=i+1;j<cntg;j++){
        if (gal.clustering[j]==1) {
          double pair_sep = sep(gal.x[i], gal.y[i], gal.z[i], gal.x[j], gal.y[j], gal.z[j]);
          if (pair_sep>sepmin && pair_sep<sepmax){
            int sum_fbr = 0;
            int sum_cov = nx*ny;
            for (int is=0;is<ny;is++) {
#pragma omp atomic
              sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is) & bw_weights(j, is)));
            }
            double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
#pragma omp atomic
            norm += gal.syst[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight;
#pragma omp atomic
            z_eff += gal.syst[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*(gal.zs[i] + gal.zs[j]) / 2.;
          }
        }
      }
    }
    if ((i+1)%10000==0) {
      cout<<"Completed "<<i<<" galaxies";
      print_time();
    }
  }
  cout<<"Finished Double Loop";
  print_time();
  
  z_eff = z_eff / norm;

  cout<<"norm: "<<norm<<endl;    
  cout<<"z_eff: "<<z_eff<<endl;

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"norm\tz_eff"<<endl<<setprecision(6)<<fixed<<norm<<"\t"<<z_eff;
  pfw.close();
  
  return (0);
}
