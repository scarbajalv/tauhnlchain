// Válido solo para ebeam = 120GeV.
// Para otras energías modificar.

// ./X idGun index

#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <iomanip>
#include <math.h>
#include <string>
#include <cmath>
#include <sstream>
#include <random>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "cpp-getEnv.h"


using namespace std;

double cdfsigmapt2(double x, double b){
	double sigmapt2return = (1./b)*log(1./(1-x));
	return sigmapt2return;
}

int main(int argc, char *argv[]){

	if (argc <  3 ){
		cout << "FATAL ERROR: ./dsprodcppvX idGun index" << endl;
		exit(EXIT_FAILURE);
	}
	
	string idGun_s_aux = argv[1];
  int idGun = stoi(idGun_s_aux);
  double m, a, b;

	// *********** AQUÍ LOS PRINCIPALES PARÁMETROS ************
	
	int nevents = 10000000;

	switch(idGun){
  	case 431:		m = 1.96834; // mass of Dmeson
								a = 16.6339541; // x-gauss
								b = 0.63424038; // pT2-exp
								break;
		case -431:	m = 1.96834; // mass of Dmeson
								a = 15.09998238; // x-gauss
								b = 0.80034524; // pT2-exp
								break;
		case 411:		m = 1.86965; // mass of Dmeson
								a = 15.1914842; // x-gauss
								b = 0.7224411; // pT2-exp
								break;
		case -411:	m = 1.86965; // mass of Dmeson
								a = 13.1260718; // x-gauss
								b = 0.92910691; // pT2-exp
								break;
		default:		cout << "FATAL ERROR: Invalid idGun" << endl;
								exit(EXIT_FAILURE);
  }

	// ********************************************************
	  
  	double ebeam = 120;

	// setup seed as time in nanosecods
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  double seed = (time_t)ts.tv_nsec;
  // set up random number function 
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  // set seed for rnd generator  
  gsl_rng_set (r, seed);

  double sigma_param = 1./sqrt(2*a);

  vector <double> row;
	vector < vector <double> > ivector;
	vector < vector <double> > ovector;	

	double mbeam = 0.968;
	double mtarget = 0.968;

	double sqrts = sqrt(pow(ebeam+mtarget,2)-(pow(ebeam,2)-pow(mbeam,2)));
	double pzcm, pt, ecm, gamma, v, pz, e, theta, phi, px, py;
	double x_aux, pT2_aux;

	// Generamos la rawdata ¡¡¡¡¡ESTO TAL VEZ NO SEA NECESARIO!!!!!
	cout<<"Generando data..."<<endl;
	for (int i=0; i<nevents; ++i){
		x_aux = gsl_ran_gaussian(r, sigma_param);
		pT2_aux = cdfsigmapt2( gsl_rng_uniform (r) , b);

		row.clear();
		pzcm = 0.5*sqrts*x_aux;
		pt = sqrt(1.*pT2_aux);
		ecm = sqrt(1.*pow(pzcm,2)+pow(pt,2)+pow(m,2));
		gamma = (pow(mtarget,2) - pow(mbeam,2) + pow(sqrts,2))/(2*sqrts*mtarget);
		v = sqrt(1-1./pow(gamma,2));
		pz = 1.*gamma*(pzcm-(-v)*ecm);
		e = sqrt(1.*pow(pz,2)+pow(pt,2)+pow(m,2));
		theta = atan(pt/pz);
		phi = gsl_rng_uniform (r)*2*3.14159265358;
		// Debug en caso de que cdfsigmapt2 se haga infinito (x=1)
		if (isnan(e)||isnan(theta)||isinf(e)||isinf(theta)){
			cout<<" nan or inf found!"<<endl;
			continue;
		}
		row.push_back(e);
		row.push_back(theta);
		row.push_back(phi);
		ovector.push_back(row);
		cout<<"\r"<<i+1;
	}	
	cout<<endl;
		
	// Exportamos dsdata-ebeam-index.dat

	string idGun_s;
  switch(idGun){
    case 431: idGun_s="431"; break;
    case -431: idGun_s="431bar"; break;
    case 411: idGun_s="411"; break;
    case -411: idGun_s="411bar"; break;
    default: break;
  }

	stringstream ofilenamess;
	ofilenamess<< getEnv("mainconfig") << "d" << idGun_s <<"data-"<<ebeam<<"-"<<argv[2]<<".dat";
	string ofilename=ofilenamess.str();

	ofstream ofile;
	ofile.open(ofilename);

	cout<<"Exportando data..."<<endl;
	for (int i=0; i<ovector.size(); ++i){
		for (int j=0; j<3; ++j){
			ofile<<setprecision(10)<<ovector[i][j]<<" ";
		}
		ofile<<endl;
		cout<<"\r"<<i+1;
	}
	cout << endl;
	cout << ofilename << endl;
	cout<<endl<<"SUCCESS!"<<endl<<endl;
	
	
	return 0;
}