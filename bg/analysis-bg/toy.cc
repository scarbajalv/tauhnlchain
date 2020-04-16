#include <iostream>
#include <fstream> 
#include <vector>

#include <stdlib.h>
#include <gsl/gsl_histogram.h>

#include <TH1.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

int main(){

	double a[2];

	size_t nbins=3;

	gsl_histogram * h[2];
	h[0] = gsl_histogram_alloc (nbins);
	h[1] = gsl_histogram_alloc (nbins);

	double xmin=0., xmax=3.;
  double xrange[nbins+1];
  for(int i=0; i < nbins+1; ++i){
    xrange[i] =  xmin + (xmax-xmin)/nbins*i;
  }

  for (int i=0; i<2; ++i){
    gsl_histogram_set_ranges (h[i], xrange, nbins+1);
    gsl_histogram_set_ranges (h[i], xrange, nbins+1);
  }

  gsl_histogram_increment(h[0],0.5);
  gsl_histogram_increment(h[0],0.5);
  gsl_histogram_increment(h[0],0.5);
  gsl_histogram_increment(h[0],1.5);
  gsl_histogram_increment(h[0],1.5);
  gsl_histogram_increment(h[0],2.5);

  for (int i=0; i<nbins; ++i){
  	cout << gsl_histogram_get(h[0],i) << endl;
  }



}