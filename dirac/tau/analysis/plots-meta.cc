// ./X idGun

// nutotal vs offaxis ( mHNL = 0.1, ... , 1.9 )

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

  // *****************************************************

  string hnlNature = "dirac"; // dirac v majorana
  string hnlMother = "tau"; // meson v tau
  
  // ******************************************************

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
  nroot_ss  << getEnv("tauhnldata") << hnlNature << "/" << hnlMother << "/" 
            << idGun_s << "/" << idGun_s << ".root";
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
  nui->SetBranchAddress("mhnl",&rmhnl);
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
  int D_mhnl = 15; // debería ser 15
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
  nEvents = 1*jEvents;

  cout << endl << "Data file:" << endl;
  cout << " "<< nroot_s << endl;
  cout << "jEvents = " << jEvents << endl;
  cout << "nEvents = " << nEvents << endl << endl;
  cout << "Realizando análisis..." << endl;
  // *********************  MAIN LOOP  **************************

  for (int i=0; i<nEvents; ++i){

		nui->GetEntry(i);
    
    if( (rdet_id==0    /*|| rdet_id==1*/   || rdet_id==2) &&
        rmhnl > 0.49 ){
      
      k = int((rmhnl-0.5)/0.1+0.0001);
      //printf("%.2f %d \n", rmhnl, k);

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

  stringstream data_ss;
  data_ss << idGun_s << "-meta.dat";
  string data_s = data_ss.str();
  char data_c[data_s.size()+1];
  strcpy(data_c, data_s.c_str());
  ofstream data(data_c);

  stringstream databar_ss;
  databar_ss << idGun_s << "-metabar.dat";
  string databar_s = databar_ss.str();
  char databar_c[databar_s.size()+1];
  strcpy(databar_c, databar_s.c_str());
  ofstream databar(databar_c);

  stringstream dataBRfactors_ss;
  dataBRfactors_ss << getEnv("tauhnlchain") << hnlNature << "/" << hnlMother 
    << "/mathbrs/brfactors" << idGun_s << ".dat"; 
  string dataBRfactors_s = dataBRfactors_ss.str();
  char dataBRfactors_c[dataBRfactors_s.size()+1];
  strcpy(dataBRfactors_c, dataBRfactors_s.c_str());
  fstream idata(dataBRfactors_s, std::ios_base::in);    
  vector <double> row;
  vector < vector <double> > v_brfactors;
  double a;
  int jcount = 0, icount = 0;
  while (idata>>a){
    jcount++;
    row.push_back(a);
    if (jcount % 2 == 0){
      v_brfactors.push_back(row);
      row.clear();
      icount++;
    }
  }

  cout << endl <<endl << "Utilizando BR factors de:" << endl 
    << " " << dataBRfactors_s << endl; 

  cout << endl << "Exportando data..." << endl;

  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        data << setprecision(16) << v_meta1[k][j][m]*v_brfactors[k][1] << " ";
      }
      data << endl;
    }
  }

  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        databar << setprecision(16) << v_meta2[k][j][m]*v_brfactors[k][1] << " ";
      }
      databar << endl;
    }
  }

  cout << data_s << endl;
  cout << databar_s << endl;
  
  cout << endl << "SUCCESS!" << endl << endl;

  return 0;
	
}

