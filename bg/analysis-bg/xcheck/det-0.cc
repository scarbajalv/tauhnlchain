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

#include "detect.h"

using namespace std;

int main(){

	double w = 7.; // width X (m)
	double h = 3.; // height Y (m)
	double l = 5.; // length Z (m)
	double z0det = 574; // z0 (m)
	double unitfactor = 1000; // 1000 for mm

	double offaxis = 21351;

	double	xidet = (-w/2+offaxis)*unitfactor, 
					xfdet = (w/2+offaxis)*unitfactor, 
					yidet = -h/2*unitfactor, 
					yfdet = h/2*unitfactor,
					zidet = z0det*unitfactor,
					zfdet = (z0det+l)*unitfactor;

	struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  double seed = (time_t)ts.tv_nsec;

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, seed);

  double xx[3];
  double x[3];
  double p[3];
  double par;
  double theta, phi;

  int counter = 0;

  int nevents = 10000000;

  for(int i=0; i<nevents; ++i){
  	xx[0] = gsl_ran_flat(r,xidet,xfdet);
  	xx[1]	=	gsl_ran_flat(r,yidet,yfdet);
  	xx[2]	=	gsl_ran_flat(r,zidet,zfdet);

  	theta =  gsl_ran_flat(r,0,3.1415926);
  	phi =  gsl_ran_flat(r,0,2*3.1415926);
  	 
  	p[0] = sin(theta)*cos(phi);
  	p[1] =	sin(theta)*sin(phi);
  	p[2] =	cos(phi);

  	par = gsl_ran_flat(r,-0,1000000);
  	x[0] = xx[0]+par*p[0];
  	x[1] = xx[1]+par*p[1];
  	x[2] = xx[2]+par*p[2];

  	if (detect(x,p,offaxis)){
  		//cout << "detected" << endl;
  		counter++;
  	}

  	/*
  	printf("%.2f %.2f %.2f  \n %.2f %.2f %.2f  \n \n"
  					,x[0],x[1],x[2],p[0],p[1],p[2]);
		*/

  }

  printf("%.4f % \n", 1.*counter/nevents*100 );
 
  return 0;


}