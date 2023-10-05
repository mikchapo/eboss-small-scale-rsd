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

  string data_input_ngc, data_input_sgc, bw_input_ngc, bw_input_sgc, dd_norm_output;
  cin>>data_input_ngc;
  cin>>data_input_sgc;
  cin>>bw_input_ngc;
  cin>>bw_input_sgc;
  cin>>dd_norm_output;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////

  cout<<"Loading data cat"<<endl;

  lrg gal_ngc;
  sfname<<data_input_ngc;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal_ngc.set_fname (fname);
  gal_ngc.fill_cat_pip();

  lrg gal_sgc;
  sfname<<data_input_sgc;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal_sgc.set_fname (fname);
  gal_sgc.fill_cat_pip();
  
  cout<<"Loading bw weights"<<endl;

  int nx = 31;
  int ny = 60;
 
  int cntg_ngc = (int)(gal_ngc.id.size());


  vector<long int> bw_ids_ngc(cntg_ngc);
  matrix<int> bw_weights_ngc(cntg_ngc, ny);
  cout<<"bw_input_ngc: "<<bw_input_ngc<<endl;

  load_bw(bw_input_ngc, ny, bw_ids_ngc, bw_weights_ngc);


  int cntg_sgc = (int)(gal_sgc.id.size());


  vector<long int> bw_ids_sgc(cntg_sgc);
  matrix<int> bw_weights_sgc(cntg_sgc, ny);
  cout<<"bw_input_sgc: "<<bw_input_sgc<<endl;

  load_bw(bw_input_sgc, ny, bw_ids_sgc, bw_weights_sgc);



  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  double D_eff_ngc = 0.;
  double D_eff_sgc = 0.;
  int i;

  cout<<"Calculating DD norm"<<endl;
  for (i=0;i<cntg_ngc;i++){
    if (gal_ngc.clustering[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
        sum_fbr+=(int)(__builtin_popcount(bw_weights_ngc(i, is)));
      }
      double iip_weight = (double)(sum_cov)/(double)(sum_fbr);
      D_eff_ngc += gal_ngc.syst[i]*gal_ngc.noz[i]*gal_ngc.fkp[i]*iip_weight;
    }
  }

  for (i=0;i<cntg_sgc;i++){
    if (gal_sgc.clustering[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
        sum_fbr+=(int)(__builtin_popcount(bw_weights_sgc(i, is)));
      }
      double iip_weight = (double)(sum_cov)/(double)(sum_fbr);
      D_eff_sgc += gal_sgc.syst[i]*gal_sgc.noz[i]*gal_sgc.fkp[i]*iip_weight;
    }
  }

  cout<<"Finished DD ";
  print_time();
  
  double DD_norm_cross = D_eff_ngc * D_eff_sgc;

  cout<<"D_eff_ngc: "<<D_eff_ngc<<" D_eff_sgc: "<<D_eff_sgc<<" DD_norm_cross: "<<DD_norm_cross<<endl;    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  sfname<<dd_norm_output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"D_eff_ngc\tD_eff_sgc\tDD_norm_cross"<<endl<<setprecision(6)<<fixed<<D_eff_ngc<<"\t"<<D_eff_sgc<<"\t"<<DD_norm_cross;
  pfw.close();
  
  return (0);
}
