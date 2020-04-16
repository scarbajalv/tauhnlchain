// ./X idGun

/* nu vs E  dado un offaxis.
 Solo cuenta eventos, No escala nada.
 nbins e_i e_f counts1 counts2 ...
            | e
  offaxis   | mu
            | tau
  (...)
*/

#include <iostream>
#include <fstream> 
#include <vector>

#include <stdlib.h>
#include <gsl/gsl_histogram.h>

#include <TH1.h>
#include <TTree.h>
#include <TFile.h>

#include "cpp-getEnv.h"

using namespace std;

int main(int argc, char *argv[]){
	
  // CONFIG NEUTRINOS
  double x0[2]={0,50};
  double x1[2]={0,50};
  double x2[2]={0,18};
  size_t nbins[3]={40,40,30};

  // CONFIG ANTINEUTRINOS
  double xbar0[2]={0,25};
  double xbar1[2]={0,45};
  double xbar2[2]={0,50};
  size_t nbinsbar[3]={40,40,40};


	// ******************************************************
	// Verificar ENERGÍA 	

  if(argc < 2){
    cout << "FATAL ERROR: ./X idGun" << endl;
    exit(EXIT_FAILURE);
  }
	
  string idGun_s_aux = argv[1];
  int idGun = stoi(idGun_s_aux);
  string idGun_s;
  switch(idGun){
    case 431: idGun_s="431"; break;
    case -431: idGun_s="431bar"; break;
    case 411: idGun_s="411"; break;
    case -411: idGun_s="411bar"; break;
    default: break;
  }

  stringstream iroot_ss;
  iroot_ss  << getEnv("hnldata") << "bg/" << idGun_s << "/" 
            << "120.root";
  string iroot_s = iroot_ss.str();
  char iroot_c[iroot_s.size()+1];
  strcpy(iroot_c, iroot_s.c_str());
  TFile f1(iroot_c, "open");  

	TTree* nui;

	f1.GetObject("nu",nui);
	
	double rebeam, rmhnl, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
  rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rmotherid;
  int mother1, auxpdg, rdet_id;

	nui->SetBranchAddress("ebeam",&rebeam);
  //nui->SetBranchAddress("mhnl",&rmhnl);
  nui->SetBranchAddress("offaxis",&roffaxis);
  nui->SetBranchAddress("det_id",&rdet_id);
  nui->SetBranchAddress("id",&rid);
  nui->SetBranchAddress("tProd",&rtProd);
  nui->SetBranchAddress("xProd",&rxProd);
  nui->SetBranchAddress("yProd",&ryProd);
  nui->SetBranchAddress("zProd",&rzProd);
  nui->SetBranchAddress("tDec",&rtDec);
  nui->SetBranchAddress("xDec",&rxDec);
  nui->SetBranchAddress("yDec",&ryDec);
  nui->SetBranchAddress("zDec",&rzDec);
  nui->SetBranchAddress("e",&re);
  nui->SetBranchAddress("px",&rpx);
  nui->SetBranchAddress("py",&rpy);
  nui->SetBranchAddress("pz",&rpz);
  nui->SetBranchAddress("pT",&rpT);
  nui->SetBranchAddress("theta",&rtheta);
  nui->SetBranchAddress("phi",&rphi);
  nui->SetBranchAddress("y",&ry);
  nui->SetBranchAddress("eta",&reta);
  //nui->SetBranchAddress("motherid",&rmotherid);
	
	double treesize = nui->GetEntries();
	
	cout << endl << "Input file: " <<  endl << " " << iroot_s << endl;
  cout << "nEvents: "<< setprecision(12) << treesize << endl;
	

  
  // toy mHNL, 41 offaxis, 3 sabores
  int D_mhnl = 1; // debería ser 15
  gsl_histogram * h[D_mhnl][41][3];
  gsl_histogram * hbar[D_mhnl][41][3];

  for (int i = 0; i < D_mhnl; ++i){ // ningún HNL
    for (int j = 0; j < 41; ++j){ // 40 offaxis
      for (int k=0; k < 3; ++k){ // 3 sabores
          h[i][j][k] = gsl_histogram_alloc (nbins[k]);
          hbar[i][j][k] = gsl_histogram_alloc (nbinsbar[k]);     
      }
    }
  }
  
  double xrange_0[nbins[0]];
  double xrange_1[nbins[1]];
  double xrange_2[nbins[2]];
  double xrangebar_0[nbinsbar[0]];
  double xrangebar_1[nbinsbar[1]];
  double xrangebar_2[nbinsbar[2]];


  
  for(int q=0; q < nbins[0]+1; ++q){
            xrange_0[q] =  x0[0] + (x0[1]-x0[0])/nbins[0]*q;
  }
  for(int q=0; q < nbins[1]+1; ++q){
            xrange_1[q] =  x1[0] + (x1[1]-x1[0])/nbins[1]*q;
  }
  for(int q=0; q < nbins[2]+1; ++q){
            xrange_2[q] =  x2[0] + (x2[1]-x2[0])/nbins[2]*q;
  }

  for(int q=0; q < nbinsbar[0]+1; ++q){
            xrangebar_0[q] =  xbar0[0] + (xbar0[1]-xbar0[0])/nbinsbar[0]*q;
  }
  for(int q=0; q < nbinsbar[1]+1; ++q){
            xrangebar_1[q] =  xbar1[0] + (xbar1[1]-xbar1[0])/nbinsbar[1]*q;
  }
  for(int q=0; q < nbinsbar[2]+1; ++q){
            xrangebar_2[q] =  xbar2[0] + (xbar2[1]-xbar2[0])/nbinsbar[2]*q;
  }

  //gsl_histogram_set_ranges_uniform (histo, XRANGE, XRANGE_Dim);

  for (int i = 0; i < D_mhnl; ++i){
    for (int j = 0; j < 41; ++j){
      for (int k=0; k < 3; ++k){
        if(k==0){
          gsl_histogram_set_ranges (h[i][j][k], xrange_0, nbins[0]+1);
          gsl_histogram_set_ranges (hbar[i][j][k], xrangebar_0, nbinsbar[0]+1);  
        }
        if(k==1){
          gsl_histogram_set_ranges (h[i][j][k], xrange_1, nbins[1]+1);
          gsl_histogram_set_ranges (hbar[i][j][k], xrangebar_1, nbinsbar[1]+1);  
        }
        if(k==2){
          gsl_histogram_set_ranges (h[i][j][k], xrange_2, nbins[2]+1);
          gsl_histogram_set_ranges (hbar[i][j][k], xrangebar_2, nbinsbar[2]+1);  
        }
      }
    }
  }

  int k = 0;
  int l = 0;

  int nevents;
  nevents = 1.*treesize;


  // *********************  MAIN LOOP  **************************

  for (int i=0; i<nevents; ++i){

		nui->GetEntry(i);

    if( ( rdet_id==0    /*|| rdet_id==1*/   || rdet_id==2 ) ){
      
      k = 0;
      
      l = int(roffaxis + 0.0001);

      //cout << k << " " << l << endl;

      switch(int(rid+0.0001)){
        case 12:
          gsl_histogram_increment(h[k][l][0],re); break;
        case 14:
          gsl_histogram_increment(h[k][l][1],re); break;
        case 16:
          gsl_histogram_increment(h[k][l][2],re); break;
        default:
          break;
      }

      switch(int(rid-0.0001)){
        case -12:
          gsl_histogram_increment(hbar[k][l][0],re); break;
        case -14:
          gsl_histogram_increment(hbar[k][l][1],re); break;
        case -16:
          gsl_histogram_increment(hbar[k][l][2],re); break;
        default:
          break;
      }

    }
		
    if( (i+1) % int(nevents/1000) == 0 || i+1 == nevents){  
      printf("%d %.1f% \r", i+1, 1.*(i+1)/nevents*100);  
    }
    fflush(stdout);
	}
  cout << endl << endl; 

  stringstream data_ss;
  stringstream databar_ss;
  data_ss << idGun_s << "-e.dat";
  databar_ss << idGun_s << "-ebar.dat";
  string data_s=data_ss.str();
  string databar_s=databar_ss.str();

  ofstream data(data_s);
  ofstream databar(databar_s); 
  
  for (int k = 0; k < D_mhnl; ++k){
    for(int l = 0; l < 41; ++l){
      for(int m = 0; m < 3; ++m){
        switch(m){
          case 0: data  << nbins[0] << " " 
                        << x0[0] << " "
                        << x0[1] << " "; break;
          case 1: data  << nbins[1] << " " 
                        << x1[0] << " "
                        << x1[1] << " "; break;
          case 2: data  << nbins[2] << " " 
                        << x2[0] << " "
                        << x2[1] << " "; break;
          default: break;
        }
        for (int i=0; i < nbins[m]; ++i){
            data << gsl_histogram_get( h[k][l][m], i)<< " ";
        }
        data << endl;
      }
    }
  }

for (int k = 0; k < D_mhnl; ++k){
    for(int l = 0; l < 41; ++l){
      for(int m = 0; m < 3; ++m){
        switch(m){
          case 0: databar   << nbinsbar[0] << " " 
                            << xbar0[0] << " "
                            << xbar0[1] << " "; break;
          case 1: databar   << nbinsbar[1] << " " 
                            << xbar1[0] << " "
                            << xbar1[1] << " "; break;
          case 2: databar   << nbinsbar[2] << " " 
                            << xbar2[0] << " "
                            << xbar2[1] << " "; break;
          default: break;
        }
        for (int i=0; i < nbinsbar[m]; ++i){
            databar << gsl_histogram_get( hbar[k][l][m], i)<< " ";
        }       
        databar << endl;
      }
    }
  }

  cout << data_s << endl;
  cout << databar_s << endl;

  cout << endl << "SUCCESS!" << endl;

  return 0;
	
}

