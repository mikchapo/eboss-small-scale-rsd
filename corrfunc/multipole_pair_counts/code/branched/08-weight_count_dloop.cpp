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
  string data_file, rand_file, output_file, fill_type;
  cin>>data_file;
  cin>>rand_file;
  cin>>output_file;
  

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  elg gal;
  sfname<<data_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal.set_fname (fname);
  
  if (fill_type.compare("weighted")==0) {
    gal.fill_cat_weighted();
  } else {
    gal.fill_cat();
  }
   
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
  double dd_pair_tot = 0, dr_pair_tot = 0, rr_pair_tot = 0;
  int i,j;

  for (i=0;i<cntg;i++){
    if (gal.fiber[i]==1) {
      for(j=i+1;j<cntg;j++) {
        if (gal.fiber[j]==1) {
          dd_pair_tot = dd_pair_tot + gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*gal.wg[i]*gal.wr[i]*gal.wz[i]*gal.syst[j]*gal.cp[j]*gal.noz[j]*gal.fkp[j]*gal.wg[j]*gal.wr[j]*gal.wz[j];
        }
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntg;
      print_time();
    }
  }

  cout<<dd_pair_tot<<endl;


  for (i=0;i<cntg;i++){
    if (gal.fiber[i]==1) {
      for(j=0;j<cntr;j++) {
        dr_pair_tot = dr_pair_tot + gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*gal.wg[i]*gal.wr[i]*gal.wz[i]*rnd.syst[j]*rnd.cp[j]*rnd.noz[j]*rnd.fkp[j];
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntg;
      print_time();
    }
  }

  cout<<dr_pair_tot<<endl;


  for (i=0;i<cntr;i++){
    for (j=i+1;j<cntr;j++){
      rr_pair_tot = rr_pair_tot + rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i]*rnd.syst[j]*rnd.cp[j]*rnd.noz[j]*rnd.fkp[j];
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntr;
      print_time();
    }
  }

  cout<<rr_pair_tot<<endl;
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output_file;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  pfw<<"DD Pairs: "<<dd_pair_tot<<endl<<"DR Pairs: "<<dr_pair_tot<<endl<<"RR Pairs: "<<rr_pair_tot<<endl;
  pfw.close();
  
  return (0);
}
