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
  string data_input, bw_file, auw_input, output;
  cin>>data_input;
  cin>>bw_file;
  cin>>auw_input;
  cin>>output;

  cout<<"Parameter Loading Complete"<<endl;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  lrg gal;
  sfname<<data_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal.set_fname (fname);
  gal.fill_cat_pip();
    
  int cntg = (int)(gal.id.size());

  int nx = 31;   // Constant value for PIP weighting
  int ny = 60;   // Number of integers used to store the bitwise weights

  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);

  load_bw(bw_file, ny, bw_ids, bw_weights);

  ang_weights auw;
  sfname<<auw_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  auw.set_fname (fname);
  auw.load_weights();

//  cout<<"AUW Loaded Values"<<endl;
//  for (int p=0; p<(int)(auw.ang.size()); p++) {
//    cout<<auw.theta[p]<<"\t"<<auw.ang[p]<<endl;
//  }


  cout<<"Catalogue Loading Complete"<<endl;

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

  double tmin = 0.01;
  double tmax = 3.;
  double st = 0.075;
  tmin = tmin * pow(10., -st/2.);
  int nBint = (int)(1./st*log10(tmax/tmin))+1;
  tmax = tmin * pow(10., nBint*st);

  cout<<"Binning Complete"<<endl;

 

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////  
  double xmin = vec_min (gal.x);
  double ymin = vec_min (gal.y);
  double zmin = vec_min (gal.z);

  cout<<"Found coord mins"<<endl;

  double xmax = vec_max (gal.x);
  double ymax = vec_max (gal.y);
  double zmax = vec_max (gal.z);

  cout<<"Found coord maxs"<<endl;

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  cout<<"Found coord ranges"<<endl;

  const int M = 300;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntg);
  matrix3d<int> label(M, M, M);
  llist (cntg, M, l, gal.x, gal.y, gal.z, lst, label);
  
  cout<<"Linked List Setup Complete"<<endl;

  

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> DD(nBinsep,nBinmu);
  int _i;
  int output_counter = 0;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (_i=0;_i<cntg;_i++){
    if (gal.clustering[_i]==1) {
      int i_1=floor(((gal.x[_i]-xmin)-sepmax)/l);
      int i_2=floor(((gal.x[_i]-xmin)+sepmax)/l);
      i_1 = max (0,i_1);
      i_2 = min (M-1,i_2);
      int j_1=floor(((gal.y[_i]-ymin)-sepmax)/l);
      int j_2=floor(((gal.y[_i]-ymin)+sepmax)/l);
      j_1 = max (0,j_1);
      j_2 = min (M-1,j_2);
      int k_1=floor(((gal.z[_i]-zmin)-sepmax)/l);
      int k_2=floor(((gal.z[_i]-zmin)+sepmax)/l);
      k_1 = max (0,k_1);
      k_2 = min (M-1,k_2);
      for (int il=i_1;il<=i_2;il++){
        for (int jl=j_1;jl<=j_2;jl++){
          for (int kl=k_1;kl<=k_2;kl++){
            int j=label(il,jl,kl);
            while (j!=0){
              if (_i<j && gal.clustering[j]==1) {
                double pair_sep = sep(gal.x[_i], gal.y[_i], gal.z[_i], gal.x[j], gal.y[j], gal.z[j]);
                if (pair_sep>sepmin && pair_sep<sepmax){
                  double r_pi = pi(gal.x[_i], gal.y[_i], gal.z[_i], gal.x[j], gal.y[j], gal.z[j]);
                  double pair_mu = r_pi / pair_sep;
		
                  if (pair_mu>mumin && pair_mu<mumax){
	            int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
	            int binmu = (int)((pair_mu-mumin)/dmu);
                  
                    double pair_sep_ang = ang_sep(gal.ra[_i], gal.dec[_i], gal.ra[j], gal.dec[j]);
                    int bint = (int)(log10(pair_sep_ang/tmin)/st);
                    double ang_weight = 1.;
                    if (bint < 0) {
                      ang_weight = auw.ang[0];
                    } else if (bint < nBint) {
                      ang_weight = auw.ang[bint];
                    }

                    int sum_fbr = 0;
                    int sum_cov = nx*ny;
                    for (int is=0;is<ny;is++) {
#pragma omp atomic
                      sum_fbr+=(int)(__builtin_popcount(bw_weights(_i,is) & bw_weights(j,is)));
                    }
                    double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
                    if (output_counter<10) {
                      cout<<"i: "<<_i<<" j: "<<j<<" ID_i: "<<gal.id[_i]<<" ID_j: "<<gal.id[j]<<" binsep: "<<binsep<<"binmu"<<binmu<<" pip_weight: "<<pip_weight<<endl;
                      cout<<"ang_weight: "<<ang_weight<<" tot_weight: "<<gal.syst[_i]*gal.noz[_i]*gal.fkp[_i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*ang_weight<<endl;
                    }
#pragma omp atomic
                    DD(binsep,binmu) = DD(binsep,binmu) + gal.syst[_i]*gal.noz[_i]*gal.fkp[_i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*ang_weight;
#pragma omp atomic
                    output_counter += 1;
                  }
                }
              }	    
              j = lst[j];
            }
          }
        }
      }
    }
    if ((_i%10000)==0){
      cout<<"I = "<<setfill('0')<<setw(6)<<_i<<"\tNg = "<<cntg;
      print_time();
    }
  }
    
  cout<<"Pair Counting Complete"<<endl;

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
      pfw<<pow(10., i*logDsep + logSepmin)<<"\t"<<j*dmu<<"\t"<<setprecision(6)<<fixed<<DD(i,j)<<endl;
    }
  }
  pfw.close();
  
  cout<<"Saving Complete"<<endl;

  return (0);
}
