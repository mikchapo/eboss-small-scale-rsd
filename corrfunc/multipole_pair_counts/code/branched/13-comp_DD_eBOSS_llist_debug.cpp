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
  string input, output, fill_type;
  cin>>input;
  cin>>output;
  cin>>fill_type;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  elg gal;
  sfname<<input;
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


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON s AND mu          //////////
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
  for (i=0;i<cntg;i++){
    cout<<i;
    if (gal.fiber[i]==1) {
      // cout<<" Position:  "<<gal.x[i]<<" "<<gal.y[i]<<" "<<gal.z[i]<<" Mins: "<<xmin<<" "<<ymin<<" "<<zmin<<" l: "<<l;
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
      // cout<<" "<<i_1<<" "<<i_2<<" "<<j_1<<" "<<j_2<<" "<<k_1<<" "<<k_2<<endl;
      for (int il=i_1;il<=i_2;il++){
        for (int jl=j_1;jl<=j_2;jl++){
          for (int kl=k_1;kl<=k_2;kl++){
            int j=label(il,jl,kl);
            while (j!=0){
              if (i<j && gal.fiber[j]==1){
                cout<<"Gal 1 Position: "<<gal.x[i]<<" "<<gal.y[i]<<" "<<gal.z[i];
                cout<<" llist indices: "<<i_1<<" "<<i_2<<" "<<j_1<<" "<<j_2<<" "<<k_1<<" "<<k_2;
                double pair_sep = sep(gal.x[i], gal.y[i], gal.z[i], gal.x[j], gal.y[j], gal.z[j]);
	        cout<<" Sep: "<<pair_sep;
                if (pair_sep>sepmin && pair_sep<sepmax){
	          double r_pi = pi(gal.x[i], gal.y[i], gal.z[i], gal.x[j], gal.y[j], gal.z[j]);
                  double pair_mu = r_pi / pair_sep;
                  cout<<" R_pi: "<<r_pi<<" Mu: "<<pair_mu<<endl;		  

                  if (pair_mu>mumin && pair_mu<mumax){
	            int binmu = (int)((pair_mu-mumin)/dmu);
                    int binsep = (int)((pair_sep-sepmin)/dsep);
#pragma omp atomic
		    DD(binsep,binmu) = DD(binsep, binmu) + gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*gal.wg[i]*gal.wr[i]*gal.wz[i]*gal.syst[j]*gal.cp[j]*gal.noz[j]*gal.fkp[j]*gal.wg[j]*gal.wr[j]*gal.wz[j];
	          }
                } else {
                  cout<<endl;
                }
	      }
	      j = lst[j];
	    }
          }
        }
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNg = "<<cntg;
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
      pfw<<i<<"\t"<<j<<"\t"<<setprecision(2)<<fixed<<DD(i,j)<<endl;
  pfw.close();
  
  return (0);
}
