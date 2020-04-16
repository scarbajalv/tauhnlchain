// ./root-merge combined-i.root

// Añade la data del archivo combined-i.root al archivo combined.root

// Agregar data de trees de un root "IN" a trees de un root "OUT".
// El root "OUT" ya debe contener trees a los cuales se agregará la data.
// Se debe conocer los nombres de los trees de interés en ambos root.

#include <iostream>
#include <fstream>


#include <TTree.h> 
#include <TFile.h>
#include <TH1.h>

#include <string>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]){
	
	// *******************************************************
	// Check ENERGY and SEED
	// OUTPUT no puede estar vacío!! (Debe ser inicializado con mergeini)

	stringstream ifile_ss;
	ifile_ss << "combined-"<<argv[1] <<".root";
	string ifile_s = ifile_ss.str();
	char ifile_c[ifile_s.size()+1];
	strcpy(ifile_c, ifile_s.c_str());
	
	ifstream ifiletest;
	ifiletest.open(ifile_s);
	//cout << ifiletest.good();
	
	if (ifiletest.good()==0){
		cout << endl << ifile_s <<" NOT FOUND"<<endl;
		exit(EXIT_FAILURE);
	}
	
	TFile f1(ifile_c, "open");

	stringstream ofile_ss;
	//ofile_ss << argv[1]  <<".root";
	ofile_ss << "combined.root";
	string ofile_s = ofile_ss.str();
	char ofile_c[ofile_s.size()+1];
	strcpy(ofile_c, ofile_s.c_str());
	
	TFile f2(ofile_c,"update"); // ¡Solo Para actualizar!
	
	// *******************************************************
	
	int id_r;
	double pz_r, pT_r;

	
	// DEFINE TREES FROM INPUT FILE:
	
	f1.cd();
	
	TTree* datai;
		
	f1.GetObject("Data",datai);

	datai->SetBranchAddress("id",&id_r);
	datai->SetBranchAddress("pz",&pz_r);
	datai->SetBranchAddress("pt",&pT_r);

		//DEFINE TREES FROM OUTPUT FILE
	
	f2.cd();
	
	TTree* datao;
		
	f2.GetObject("Data",datao);

	datao->SetBranchAddress("id",&id_r);
	datao->SetBranchAddress("pz",&pz_r);
	datao->SetBranchAddress("pt",&pT_r);
	
			
	// UPDATE TREES
	
	for (int i=0; i<datai->GetEntries();++i){
		datai->GetEntry(i); //almacenar e(i) en var1
		datao->Fill(); // almacenar var1 en tDp
	}

	
		
	//Write trees (removing history).
	f2.cd();
	
	datao->Write("",2);

	return 0;


	
}
