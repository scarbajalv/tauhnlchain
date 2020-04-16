// ./X idGun

// nutotal vs offaxis 
// Solo cuenta eventos, no escala nada.

#include <iostream> 
#include <fstream> 
#include <vector>
#include <string>
#include <sstream>

#include <stdlib.h>
#include <gsl/gsl_histogram.h>

#include <TH1.h>
#include <TTree.h>
#include <TFile.h>

#include "cpp-getEnv.h"

using namespace std;

int main(int argc, char *argv[]){
	
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
    case 421: idGun_s="421"; break;
    case -421: idGun_s="421bar"; break;
    default: break;
  }

	stringstream nroot_ss;
  nroot_ss  << getEnv("hnldata") << "bg/" 
            << idGun_s << "/120.root";
  string nroot_s = nroot_ss.str();
  char nroot_c[nroot_s.size()+1];
  strcpy(nroot_c, nroot_s.c_str());
  TFile f1(nroot_c, "open");

	TTree* nui;

	f1.GetObject("nu",nui);
	
	double rebeam, rmhnl, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
  rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex;
  int auxpdg, rdet_id, rmotherid;

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
	
	int jEvents = nui->GetEntries();
	
	// 0-mHNL y 41 offaxis [0, ..., 40]
  int D_mhnl = 1; // Dimension de mHNL 
  double v_meta1[D_mhnl][41][3];
  double v_meta2[D_mhnl][41][3];


  for( int i=0; i < D_mhnl; ++i){
    for (int j=0; j < 41; ++j){
      for(int k=0; k < 3; ++k){
        v_meta1[i][j][k]=0;
        v_meta2[i][j][k]=0;
      }
    }
  }

  int k = 0;
  int l = 0;

  int nEvents;
  nEvents = jEvents;

  cout << endl << "Data file:" << endl;
  cout << " "<< nroot_s << endl;
  cout << "jEvents = " << jEvents << endl;
  cout << "nEvents = " << nEvents << endl << endl;
  cout << "Realizando análisis..." << endl;
  // *********************  MAIN LOOP  **************************

  for (int i=0; i<nEvents; ++i){

		nui->GetEntry(i);
    
    if( rdet_id==0    || rdet_id==1   || rdet_id==2 ){
      
      k = 0;
      l = int(roffaxis+0.000001);

      switch(int(rid+0.0001)){
        case 12:
          v_meta1[k][l][0]++; 
          break;
        case 14:
          v_meta1[k][l][1]++; 
          break;
        case 16:
          v_meta1[k][l][2]++; 
          break;
        default:
          break;
      }

      switch(int(rid-0.0001)){
        case -12:
          v_meta2[k][l][0]++; 
          break;
        case -14:
          v_meta2[k][l][1]++; 
          break;
        case -16:
          v_meta2[k][l][2]++; 
          break;
        default:
          break;
      }     
            
    }
		
    if( (i+1) % int(nEvents/1000) == 0 || i+1 == nEvents){  
      printf("%d %.1f% \r", i+1, 1.*(i+1)/nEvents*100);  
    }
    fflush(stdout);
	}
  cout << endl << endl;

  stringstream data_ss;
  data_ss << idGun_s << "-meta.dat";
  string data_s = data_ss.str();
  ofstream data(data_s);

  stringstream databar_ss;
  databar_ss << idGun_s << "-metabar.dat";
  string databar_s = databar_ss.str();
  ofstream databar(databar_s);
  
  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        data << v_meta1[k][j][m] << " ";
      }
      data << endl;
    }
  }

  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        databar << v_meta2[k][j][m] << " ";
      }
      databar << endl;
    }
  }

  cout << data_s << endl;
  cout << databar_s << endl;
  
  cout << endl << "SUCCESS!" << endl << endl;

  return 0;
	
}

