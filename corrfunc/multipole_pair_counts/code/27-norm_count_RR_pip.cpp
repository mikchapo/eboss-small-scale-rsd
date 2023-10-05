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
  cout<<"Reading Input"<<endl;

  string rand_file, output_file;
  cin>>rand_file;
  cin>>output_file;


  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  cout<<"Loading rand cat"<<endl;

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
  double R_eff_1 = 0.;
  double R_eff_2 = 0.;
  int i;
  
  cout<<"Calculating RR norm"<<endl;
  for (i=0;i<cntr;i++){
    R_eff_1 += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
    R_eff_2 += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i]*rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
  }

  double RR_norm = (R_eff_1*R_eff_1 - R_eff_2) / 2.;

  cout<<"R_eff_1: "<<R_eff_1<<" R_eff_2: "<<R_eff_2<<" RR_norm: "<<RR_norm<<endl;
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"R_eff_1\tR_eff_2\tRR_norm"<<endl<<setprecision(6)<<fixed<<R_eff_1<<"\t"<<R_eff_2<<"\t"<<RR_norm;
  pfw.close();
  
  return (0);
}
