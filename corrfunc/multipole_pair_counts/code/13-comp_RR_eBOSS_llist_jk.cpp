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
  
  int nReg;
  cin>>nReg;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rand;
  sfname<<input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  rand.set_fname (fname);
  rand.fill_cat_jk();
    
  int cntr = (int)(rand.ra.size());


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double sepmin = 0.1;
  double sepmax = 60.;
  double logSepmin = log10(sepmin);
  double logSepmax = log10(sepmax);
  int nBinsep = 9;
  double dsep = ((sepmax-sepmin)/(double)nBinsep);
  double logDsep = ((logSepmax - logSepmin)/(double)nBinsep);
 
  double mumin = 0.;
  double mumax = 1.;
  double dmu = 0.01;
  int nBinmu = (int)((mumax-mumin)/dmu);
 

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (rand.x);
  double ymin = vec_min (rand.y);
  double zmin = vec_min (rand.z);

  double xmax = vec_max (rand.x);
  double ymax = vec_max (rand.y);
  double zmax = vec_max (rand.z);

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  const int M = 300;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntr);
  matrix3d<int> label(M, M, M);
  llist (cntr, M, l, rand.x, rand.y, rand.z, lst, label);
    

   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix4d<double> RR(nBinsep, nBinmu, nReg, nReg);
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
	          int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
	          int binmu = (int)((pair_mu-mumin)/dmu);
                  int reg1, reg2;
                  if (rand.jk[i]==-1 && rand.jk[j]==-1) {
                    reg1 = nReg-1;
                    reg2 = nReg-1;
                  } else if (rand.jk[i]==-1) {
                    reg1 = nReg-1;
                    reg2 = rand.jk[j];
                  } else if (rand.jk[j]==-1) {
                    reg1 = rand.jk[i];
                    reg2 = nReg-1;
                  } else {
                    reg1 = rand.jk[i];
                    reg2 = rand.jk[j];
                  }

//		  if (binsep==0 && binmu==1 && reg1==137 && reg2==137) {
//                    cout<<"Bin Match:"<<endl;
//                    cout<<rand.ra[i]<<"\t"<<rand.dec[i]<<"\t"<<rand.x[i]<<"\t"<<rand.y[i]<<"\t"<<rand.z[i]<<"\t"<<rand.syst[i]<<"\t"<<rand.cp[i]<<"\t"<<rand.noz[i]<<"\t"<<rand.fkp[i]<<"\t"<<rand.jk[i]<<endl;
//                    cout<<rand.ra[j]<<"\t"<<rand.dec[j]<<"\t"<<rand.x[j]<<"\t"<<rand.y[j]<<"\t"<<rand.z[j]<<"\t"<<rand.syst[j]<<"\t"<<rand.cp[j]<<"\t"<<rand.noz[j]<<"\t"<<rand.fkp[j]<<"\t"<<rand.jk[j]<<endl;
//                    cout<<pair_sep<<"\t"<<r_pi<<"\t"<<pair_mu<<"\t"<<endl;
//                  }

//		  if (abs(rand.syst[i]*rand.cp[i]*rand.noz[i]*rand.fkp[i]*rand.syst[j]*rand.cp[j]*rand.noz[j]*rand.fkp[j]-0.267519)<0.000001) {
//                    cout<<"Weight Match:"<<endl;
//                    cout<<rand.ra[i]<<"\t"<<rand.dec[i]<<"\t"<<rand.x[i]<<"\t"<<rand.y[i]<<"\t"<<rand.z[i]<<"\t"<<rand.syst[i]<<"\t"<<rand.cp[i]<<"\t"<<rand.noz[i]<<"\t"<<rand.fkp[i]<<"\t"<<rand.jk[i]<<endl;
//                    cout<<rand.ra[j]<<"\t"<<rand.dec[j]<<"\t"<<rand.x[j]<<"\t"<<rand.y[j]<<"\t"<<rand.z[j]<<"\t"<<rand.syst[j]<<"\t"<<rand.cp[j]<<"\t"<<rand.noz[j]<<"\t"<<rand.fkp[j]<<"\t"<<rand.jk[j]<<endl;
//                    cout<<pair_sep<<"\t"<<r_pi<<"\t"<<pair_mu<<"\t"<<binsep<<"\t"<<binmu<<"\t"<<reg1<<"\t"<<reg2<<endl;
//                  }

#pragma omp atomic
	          RR(binsep, binmu, reg1, reg2) = RR(binsep, binmu, reg1, reg2) + rand.syst[i]*rand.cp[i]*rand.noz[i]*rand.fkp[i]*rand.syst[j]*rand.cp[j]*rand.noz[j]*rand.fkp[j];
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
  for (int i=0;i<nBinsep;i++) {
    for (int j=0;j<nBinmu;j++) {
      for (int k=0;k<nReg;k++) {
        for (int l=0;l<nReg;l++) {
          if (RR(i,j,k,l)!=0) {
            pfw<<pow(10., i*logDsep + logSepmin)<<"\t"<<j*dmu<<"\t"<<k<<"\t"<<l<<"\t"<<setprecision(6)<<fixed<<RR(i,j,k,l)<<endl;
          }
        }
      }
    }
  }
  pfw.close();
  
  return (0);
}
