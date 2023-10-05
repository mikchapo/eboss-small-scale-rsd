#include "catalogue.h"








void lrg::fill_cat_basic(){

  ifstream pf;
  string temp;  
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
	if (temp[0]!='#' && temp[0]!='r'){
	  vector<double> vec;
	  char *p;
	  char* temp1 = strdup(temp.c_str());
	  p = strtok (temp1," W P ,\t\n");
	  while (p!= NULL){
	    vec.push_back(atof(p));
	    p = strtok (NULL," W P ,\t");
	  }
	  id.push_back((long int)(vec[0]));
	  ra.push_back(vec[1]);
	  dec.push_back(vec[2]);
	  zs.push_back(vec[3]);
	  syst.push_back(vec[4]);
	  cp.push_back(vec[5]);
	  noz.push_back(vec[6]);
	  fkp.push_back(vec[7]);
	  x.push_back(0.);
	  y.push_back(0.);
	  z.push_back(0.);
          jk.push_back(0.);
          fiber.push_back(1);
          clustering.push_back(1);
	  vec.clear();
	  free (p);
	  free (temp1);
	}
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void lrg::fill_cat(){

  ifstream pf;
  string temp;
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          id.push_back((long int)(vec[0]));
          ra.push_back(vec[1]);
          dec.push_back(vec[2]);
          zs.push_back(vec[3]);
	  x.push_back(vec[4]);
	  y.push_back(vec[5]);
	  z.push_back(vec[6]);
          syst.push_back(vec[7]);
          cp.push_back(vec[8]);
          noz.push_back(vec[9]);
          fkp.push_back(vec[10]);
          jk.push_back(0.);
          fiber.push_back(1);
          clustering.push_back(1);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void lrg::fill_cat_jk(){

  ifstream pf;
  string temp;
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          id.push_back((long int)(vec[0]));
          ra.push_back(vec[1]);
          dec.push_back(vec[2]);
          zs.push_back(vec[3]);
	  x.push_back(vec[4]);
	  y.push_back(vec[5]);
	  z.push_back(vec[6]);
          syst.push_back(vec[7]);
          cp.push_back(vec[8]);
          noz.push_back(vec[9]);
          fkp.push_back(vec[10]);
          jk.push_back(vec[11]);
          fiber.push_back(1);
          clustering.push_back(1);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}



void lrg::fill_cat_pip(){

  ifstream pf;
  string temp;
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          id.push_back((long int)(vec[0]));
          ra.push_back(vec[1]);
          dec.push_back(vec[2]);
          zs.push_back(vec[3]);
	  x.push_back(vec[4]);
	  y.push_back(vec[5]);
	  z.push_back(vec[6]);
          syst.push_back(vec[7]);
          cp.push_back(vec[8]);
          noz.push_back(vec[9]);
          fkp.push_back(vec[10]);
          jk.push_back(0.);
          fiber.push_back((int)(vec[11]));
          clustering.push_back((int)(vec[12]));
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void randcat::fill_cat_basic(){

  ifstream pf;
  
  
  string temp;
  
  
  
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
	if (temp[0]!='#' && temp[0]!='r'){
	  vector<double> vec;
	  char *p;
	  char* temp1 = strdup(temp.c_str());
	  p = strtok (temp1," W P ,\t\n");
	  while (p!= NULL){
	    vec.push_back(atof(p));
	    p = strtok (NULL," W P ,\t");
	  }
	  ra.push_back(vec[0]);
	  dec.push_back(vec[1]);
          syst.push_back(vec[2]);
          cp.push_back(vec[3]);
          noz.push_back(vec[4]);
          fkp.push_back(vec[5]);
	  zs.push_back(0.);
	  x.push_back(0.);
          y.push_back(0.);
	  z.push_back(0.);
          jk.push_back(0.);
	  vec.clear();
	  free (p);
	  free (temp1);
	}
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void randcat::fill_cat(){

  ifstream pf;


  string temp;



  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          ra.push_back(vec[0]);
          dec.push_back(vec[1]);
          zs.push_back(vec[2]);
          x.push_back(vec[3]);
          y.push_back(vec[4]);
          z.push_back(vec[5]);
          syst.push_back(vec[6]);
          cp.push_back(vec[7]);
          noz.push_back(vec[8]);
          fkp.push_back(vec[9]);
          jk.push_back(0.);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void randcat::fill_cat_jk(){

  ifstream pf;


  string temp;



  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          ra.push_back(vec[0]);
          dec.push_back(vec[1]);
          zs.push_back(vec[2]);
          x.push_back(vec[3]);
          y.push_back(vec[4]);
          z.push_back(vec[5]);
          syst.push_back(vec[6]);
          cp.push_back(vec[7]);
          noz.push_back(vec[8]);
          fkp.push_back(vec[9]);
          jk.push_back(vec[10]);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void ang_weights::load_weights_basic(){

  ifstream pf;
  
  
  string temp;
  
  
  
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
	if (temp[0]!='#' && temp[0]!='r'){
	  vector<double> vec;
	  char *p;
	  char* temp1 = strdup(temp.c_str());
	  p = strtok (temp1," W P ,\t\n");
	  while (p!= NULL){
	    vec.push_back(atof(p));
	    p = strtok (NULL," W P ,\t");
	  }
	  theta.push_back(vec[0]);
	  ang.push_back(1.);
	  vec.clear();
	  free (p);
	  free (temp1);
	}
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void ang_weights::load_weights(){

  ifstream pf;
  
  
  string temp;
  
  
  
  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
	if (temp[0]!='#' && temp[0]!='r'){
	  vector<double> vec;
	  char *p;
	  char* temp1 = strdup(temp.c_str());
	  p = strtok (temp1," W P ,\t\n");
	  while (p!= NULL){
	    vec.push_back(atof(p));
	    p = strtok (NULL," W P ,\t");
	  }
	  theta.push_back(vec[0]);
	  ang.push_back(vec[1]);
	  vec.clear();
	  free (p);
	  free (temp1);
	}
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}



