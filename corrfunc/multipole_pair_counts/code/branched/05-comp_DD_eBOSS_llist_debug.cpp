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
  string data_input, output;
  cin>>data_input;
  cin>>output;


  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  lrg gal;
  sfname<<data_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal.set_fname (fname);
  gal.fill_cat();
    
  int cntg = (int)(gal.id.size());


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
  vector<int> lst(cntg);
  matrix3d<int> label(M, M, M);
  llist (cntg, M, l, gal.x, gal.y, gal.z, lst, label);
    
  double xmin = vec_min (gal.x);
  double ymin = vec_min (gal.y);
  double zmin = vec_min (gal.z);

   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> DD(nBinsep,nBinmu);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=299;i<1000;i++){
    cout<<"Start Loop"<<endl;
    int i_1=floor(((gal.x[i]-xmin)-sepmax)/l);
    int i_2=floor(((gal.x[i]-xmin)+sepmax)/l);
    i_1 = max (0,i_1);
    i_2 = min (M-1,i_2);
    int j_1=floor(((gal.y[i]-ymin)-sepmax)/l);
    int j_2=floor(((gal.y[i]-ymin)+sepmax)/l);
    j_1 = max (0,j_1);
    j_2 = min (M-1,j_2);
    int k_1=floor(((gal.z[i]-zmin)-sepmax)/l);
    int k_2=floor(((gal.z[i]-zmin)+sepmax)/l);
    k_1 = max (0,k_1);
    k_2 = min (M-1,k_2);
    cout<<"Calc indices"<<endl;
    for (int il=i_1;il<=i_2;il++){
      for (int jl=j_1;jl<=j_2;jl++){
        for (int kl=k_1;kl<=k_2;kl++){
          int j=label(il,jl,kl);
          // cout<<"j: "<<j<<" il: "<<il<<" jl: "<<" kl: "<<kl<<endl;
          while (j!=0){
            if (i<j) {
              cout<<"i<j"<<endl;
              double pair_sep = sep(gal.x[i], gal.y[i], gal.z[i], gal.x[j], gal.y[j], gal.z[j]);
              cout<<"Pair Sep: "<<pair_sep<<endl;
              if (pair_sep>sepmin && pair_sep<sepmax){
                cout<<"Correct Sep"<<endl;
                double r_pi = pi(gal.x[i], gal.y[i], gal.z[i], gal.x[j], gal.y[j], gal.z[j]);
                double pair_mu = r_pi / pair_sep;
		cout<<"r_pi: "<<r_pi<<" mu: "<<pair_mu<<endl;
                if (pair_mu>mumin && pair_mu<mumax){
                  cout<<"Correct mu"<<endl;
	          int binsep = (int)((pair_sep-sepmin)/dsep);
	          int binmu = (int)((pair_mu-mumin)/dmu);
#pragma omp atomic
	          DD(binsep,binmu) = DD(binsep,binmu) + gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.cp[j]*gal.noz[j]*gal.fkp[j];
                }
              }
	    j = lst[j];
	    }
          }
        }
      }
    }
    if ((i%1)==0){
      cout<<"I = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntg;
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
      pfw<<i*dsep<<"\t"<<j*dmu<<"\t"<<setprecision(2)<<fixed<<DD(i,j)<<endl;
  pfw.close();
  
  return (0);
}
