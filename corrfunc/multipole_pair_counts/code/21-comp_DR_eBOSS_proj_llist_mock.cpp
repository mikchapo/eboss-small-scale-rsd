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


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double r_perp_min = 0.1;
  double r_perp_max = 60.;
  double log_r_perp_min = log10(r_perp_min);
  double log_r_perp_max = log10(r_perp_max);
  int n_bin_perp = 180;
  double dr_perp = ((r_perp_max-r_perp_min)/(double)n_bin_perp);
  double log_dr_perp = ((log_r_perp_max - log_r_perp_min)/(double)n_bin_perp);
 
  double r_pi_min = 0.;
  double r_pi_max = 85.;
  int n_bin_pi = 85;
  double dr_pi = (r_pi_max - r_pi_min) / (double)n_bin_pi;


 

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (rnd.x);
  double ymin = vec_min (rnd.y);
  double zmin = vec_min (rnd.z);

  double xmax = vec_max (rnd.x);
  double ymax = vec_max (rnd.y);
  double zmax = vec_max (rnd.z);

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  const int M = 300;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntr);
  matrix3d<int> label(M, M, M);
  llist (cntr, M, l, rnd.x, rnd.y, rnd.z, lst, label);

  double rmax = sqrt(r_perp_max*r_perp_max + r_pi_max*r_pi_max) + 10.;
   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> DR(n_bin_perp, n_bin_pi);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    if (gal.zs[i]>=0.6 && gal.zs[i]<=1.0) {
      int i_1=floor(((gal.x[i]-xmin)-rmax)/l);
      int i_2=floor(((gal.x[i]-xmin)+rmax)/l);
      i_1 = max (0,i_1);
      i_2 = min (M-1,i_2);
      int j_1=floor(((gal.y[i]-ymin)-rmax)/l);
      int j_2=floor(((gal.y[i]-ymin)+rmax)/l);
      j_1 = max (0,j_1);
      j_2 = min (M-1,j_2);
      int k_1=floor(((gal.z[i]-zmin)-rmax)/l);
      int k_2=floor(((gal.z[i]-zmin)+rmax)/l);
      k_1 = max (0,k_1);
      k_2 = min (M-1,k_2);
      for (int il=i_1;il<=i_2;il++){
        for (int jl=j_1;jl<=j_2;jl++){
          for (int kl=k_1;kl<=k_2;kl++){
            int j=label(il,jl,kl);
            while (j!=0){
              if (rnd.zs[j]>=0.6 && rnd.zs[j]<=1.0) {
                double r_pi = pi(gal.x[i], gal.y[i], gal.z[i], rnd.x[j], rnd.y[j], rnd.z[j]);
                if (r_pi>r_pi_min && r_pi<r_pi_max){
                  double r_perp = rp(gal.x[i], gal.y[i], gal.z[i], rnd.x[j], rnd.y[j], rnd.z[j], r_pi);
		
                  if (r_perp>r_perp_min && r_perp<r_perp_max){
	            int bin_perp = (int)((log10(r_perp)-log_r_perp_min)/log_dr_perp);
	            int bin_pi = (int)((r_pi-r_pi_min)/dr_pi);
#pragma omp atomic
	            DR(bin_perp, bin_pi) = DR(bin_perp, bin_pi) + gal.fkp[i]*rnd.fkp[j];
	          }
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
  for (int i=0;i<n_bin_perp;i++)
    for (int j=0;j<n_bin_pi;j++)
      pfw<<pow(10., i*log_dr_perp + log_r_perp_min)<<"\t"<<j*dr_pi<<"\t"<<setprecision(6)<<fixed<<DR(i,j)<<endl;
  pfw.close();
  
  return (0);
}
