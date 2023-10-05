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
  string input, output;
  cin>>input;
  cin>>output;
  

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rand;
  sfname<<input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  rand.set_fname (fname);
  rand.fill_cat();
    
  int cntr = (int)(rand.ra.size());


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double sepmin = 0.;
  double sepmax = 250.;
  double dsep = 1.;
  int nBinsep = (int)((sepmax-sepmin)/dsep);
 
  double mumin = 0.;
  double mumax = 1.;
  double dmu = 0.01;
  int nBinmu = (int)((mumax-mumin)/dmu);


 

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  const int M = 200;
  const double l = 2048./((double)(M));
  vector<int> lst(cntr);
  matrix3d<int> label(M, M, M);
  llist (cntr, M, l, rand.x, rand.y, rand.z, lst, label);
    
  double xmin = vec_min (rand.x);
  double ymin = vec_min (rand.y);
  double zmin = vec_min (rand.z);

   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> RR(nBinsep,nBinmu);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntr;i++){
    int i_1=floor(((rand.x[i]-xmin)-sepmax)/l);
    int i_2=floor(((rand.x[i]-xmin)+sepmax)/l);
    i_1 = max (0,i_1);
    i_2 = min (M-1,i_2);
    int j_1=floor(((rand.y[i]-ymin)-sepmax)/l);
    int j_2=floor(((rand.y[i]-ymin)+sepmax)/l);
    j_1 = max (0,j_1);
    j_2 = min (M-1,j_2);
    int k_1=floor(((rand.z[i]-zmin)-sepmax)/l);
    int k_2=floor(((rand.z[i]-zmin)+sepmax)/l);
    k_1 = max (0,k_1);
    k_2 = min (M-1,k_2);
    for (int il=i_1;il<=i_2;il++){
      for (int jl=j_1;jl<=j_2;jl++){
        for (int kl=k_1;kl<=k_2;kl++){
          int j=label(il,jl,kl);
          while (j!=0){
            if (i<j){
              double pair_sep = sep(rand.x[i], rand.y[i], rand.z[i], rand.x[j], rand.y[j], rand.z[j]);
              if (pair_sep>sepmin && pair_sep<sepmax){
	        double r_pi = pi(rand.x[i], rand.y[i], rand.z[i], rand.x[j], rand.y[j], rand.z[j]);
                double pair_mu = r_pi / pair_sep;
		  
	        if (pair_mu>mumin && pair_mu<mumax){
	          int binsep = (int)((pair_sep-sepmin)/dsep);
	          int binmu = (int)((pair_mu-mumin)/dmu);
#pragma omp atomic
	          RR(binsep,binmu) = RR(binsep,binmu) + rand.syst[i]*rand.cp[i]*rand.noz[i]*rand.fkp[i]*rand.syst[j]*rand.cp[j]*rand.noz[j]*rand.fkp[j];
	        }
              }
	    }
	  j = lst[j];
	  }
        }
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNr = "<<cntr;
      print_time();
    }
  }
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  pfw.open (fname.c_str(), ios::out);
  for (int i=0;i<nBinsep;i++)
    for (int j=0;j<nBinmu;j++)
      pfw<<i*dsep<<"\t"<<j*dmu<<"\t"<<setprecision(2)<<fixed<<RR(i,j)<<endl;
  pfw.close();
  
  return (0);
}
