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
  string data_input, rand_input, output;
  cin>>data_input;
  cin>>rand_input;
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

  cout<<"Data Check:"<<endl;
  cout<<"LRD_ID\tRA\t\tDEC\t\tZS\t\tX\t\tY\t\tZ\t\tSYST\t\tCP\t\tNOZ\t\tFKP"<<endl;
  for (int i=0; i<20; i++) {
    cout<<gal.id[i]<<"\t"<<gal.ra[i]<<"\t\t"<<gal.dec[i]<<"\t\t"<<gal.zs[i]<<"\t\t"<<gal.x[i]<<"\t\t"<<gal.y[i]<<"\t\t"<<gal.z[i]<<"\t\t"<<gal.syst[i]<<"\t\t"<<gal.cp[i]<<"\t\t"<<gal.noz[i]<<"\t\t"<<gal.fkp[i]<<endl;
  }


  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rnd;
  sfname<<rand_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  rnd.set_fname (fname);
  rnd.fill_cat();
    
  int cntr = (int)(rnd.ra.size());


  cout<<"Rand Check:"<<endl;
  cout<<"RA\t\tDEC\t\tZS\t\tX\t\tY\t\tZ\t\tSYST\t\tCP\t\tNOZ\t\tFKP"<<endl;
  for (int i=0; i<20; i++) {
    cout<<rnd.ra[i]<<"\t\t"<<rnd.dec[i]<<"\t\t"<<rnd.zs[i]<<"\t\t"<<rnd.x[i]<<"\t\t"<<rnd.y[i]<<"\t\t"<<rnd.z[i]<<"\t\t"<<rnd.syst[i]<<"\t\t"<<rnd.cp[i]<<"\t\t"<<rnd.noz[i]<<"\t\t"<<rnd.fkp[i]<<endl;
  }

  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double sepmin = 0.1;
  double sepmax = 60.;
  double logSepmin = log10(sepmin);
  double logSepmax = log10(sepmax);
  int nBinsep = 180;
  double dsep = ((sepmax-sepmin)/(double)nBinsep);
  double logDsep = ((logSepmax - logSepmin)/(double)nBinsep);
 
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
  llist (cntr, M, l, rnd.x, rnd.y, rnd.z, lst, label);
    
  double xmin = vec_min (rnd.x);
  double ymin = vec_min (rnd.y);
  double zmin = vec_min (rnd.z);

   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> DR(nBinsep,nBinmu);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
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
    for (int il=i_1;il<=i_2;il++){
      for (int jl=j_1;jl<=j_2;jl++){
        for (int kl=k_1;kl<=k_2;kl++){
          int j=label(il,jl,kl);
          while (j!=0){
            double pair_sep = sep(gal.x[i], gal.y[i], gal.z[i], rnd.x[j], rnd.y[j], rnd.z[j]);
            if (pair_sep>sepmin && pair_sep<sepmax){
              double r_pi = pi(gal.x[i], gal.y[i], gal.z[i], rnd.x[j], rnd.y[j], rnd.z[j]);
              double pair_mu = r_pi / pair_sep;
		
              if (pair_mu>mumin && pair_mu<mumax){
	        int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
	        int binmu = (int)((pair_mu-mumin)/dmu);
#pragma omp atomic
	        DR(binsep,binmu) = DR(binsep,binmu) + gal.syst[i]*gal.cp[i]*gal.noz[i]*gal.fkp[i]*rnd.syst[j]*rnd.cp[j]*rnd.noz[j]*rnd.fkp[j];
	      }
            }
	    j = lst[j];
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
      pfw<<pow(10., i*logDsep + logSepmin)<<"\t"<<j*dmu<<"\t"<<setprecision(6)<<fixed<<DR(i,j)<<endl;
  pfw.close();
  
  return (0);
}
