#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <algorithm>

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
  string data_file, rand_file, output_file;
  cin>>data_file;
  cin>>rand_file;
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
  gal.fill_cat();
  
 
  int cntg = (int)(gal.id.size());


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
  double D_eff_1 = 0.;
  double D_eff_2 = 0.;
  double R_eff_1 = 0.;
  double R_eff_2 = 0.;
  int i,j;

  for (i=0;i<cntg;i++){
    D_eff_1 += gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i];
    D_eff_2 += gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i];
  }

  for (i=0;i<cntr;i++){
    R_eff_1 += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
    R_eff_2 += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i]*rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
  }

  double DD_norm = (D_eff_1*D_eff_1 - D_eff_2) / 2.;
  double DR_norm = (D_eff_1*R_eff_1);
  double RR_norm = (R_eff_1*R_eff_1 - R_eff_2) / 2.;

  cout<<"DD_norm: "<<DD_norm<<endl;
  cout<<"DR_norm: "<<DR_norm<<endl;
  cout<<"RR_norm: "<<RR_norm<<endl;
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"DD_norm: "<<setprecision(6)<<fixed<<DD_norm<<endl<<"DR_norm: "<<DR_norm<<endl<<"RR_norm: "<<RR_norm<<endl;
  pfw.close();
  
  return (0);
}
