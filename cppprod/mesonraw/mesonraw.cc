// ./X nEvents(M) seed
// 	nEvents(M): Número de POTs en millones.
// 	seed: RandomSeed para PYTHIA.
// Output: combined-seed.root
#include <iostream>
#include <fstream>
#include <chrono>

#include<TTree.h>
#include<TFile.h>

#include<Pythia8/Pythia.h>

using namespace Pythia8;


int main(int argc, char* argv[]){

	// Cerrar en caso falten argumentos.
	if(argc < 1+2){
		cout << "./X nEvents(M) seed" << endl;
		exit(EXIT_FAILURE);
	}

	// *********** PARÁMETROS IMPORTANTES *************

	// Tabla de PDGcodes de partículas a almacenar
	// (sus antipartículas se almacenan automáticamente)
	int idtable[]={411,431,421,15};

	// ************************************************

	//Definir variables a partir de argumentos
	string seed_s = argv[2];
	int seed = stoi(seed_s);
	
	string nEvents_s = argv[1];
	int nEvents = stod(nEvents_s)*1e6;
	
	// **** ROOT CONFIG ****

	// String para nombre del archivo
	stringstream file_ss;
	file_ss << "combined-" << seed << ".root";
	string file_s = file_ss.str();
	char file_c[file_s.size()+1];
	strcpy(file_c, file_s.c_str());

	TFile f1(file_c,"recreate");
	
	int id_r;
	double pz_r, pt_r;
	
	TTree Data("Data","Data");

	Data.Branch("id",&id_r,"id_r/I");
	Data.Branch("pz",&pz_r,"pz_r/D");
	Data.Branch("pt",&pt_r,"pt_r/D");

	// **** PYTHIA CONFIG ****

	Pythia pythia;
	pythia.readString("Random:setSeed = on");

	// String para seed de Pythia
	stringstream seedconfig_ss;
	seedconfig_ss << "Random:seed = "<< seed;
	string seedconfig_s = seedconfig_ss.str();

	pythia.readString(seedconfig_s);
	pythia.readString("SoftQCD:nonDiffractive = on");
	pythia.readString("Beams:frameType = 1");
	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");	
	pythia.readString("Beams:eCM = 15.3034");

	pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberCount = 0");

	pythia.init();

	// Inicializar contadores de tiempo
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();

	// **** PYTHIA Main Event Loop ****

	for (int iEvent = 0; iEvent < nEvents; ++iEvent){ // Main Loop Begins

		// Imprimir configuración (DEBUG)
		if(iEvent == 0){
			printf("MAIN PARAMETERS:\n");
			printf("  nEvents = %d\n",nEvents);
			printf("  seed = %d\n\n",seed);

		}

		// Generar el evento. Si PYTHIA colapsa, interrumpir el Main Loop.
		if(!pythia.next()){
			cout << "FATAL ERROR: Event generator failed." << endl;
			break;
		}

		// "Particle Loop" sobre todas las partículas del evento.
		for (int i=0; i<pythia.event.size(); ++i){ // Particle Loop Begins

			bool idcheck=false;
			
			// Verificar si la partícula nos interesa.
			for ( int j=0; j<sizeof(idtable)/sizeof(idtable[0]); ++j ){
				if (		pythia.event[i].id() == idtable[j] 
						|| -pythia.event[i].id() == idtable[j]  ){
					idcheck = true;
				}

			}
			// Si la partícula nos interesa, guardar su data en el Tree.
			if (idcheck){
					id_r = pythia.event[i].id();
					pz_r = pythia.event[i].pz();
					pt_r = pythia.event[i].pT();
					Data.Fill();
			}
			
		}// Particle Loop Ends

		// Imprimir el porcentaje completado y un estimado del tiempo restante.
		if(iEvent % (nEvents/1000) == 0){
			end = std::chrono::steady_clock::now();
			int sleft = int(1.*std::chrono::duration_cast<std::chrono::seconds>
				(end - start).count()/iEvent*(nEvents-iEvent));				
			cout << "\r" << fixed  <<setprecision(1) << 1.*iEvent/nEvents*100 
						<< "% / Tiempo restante ~ " 
						<< sleft/3600 << "h" 
						<< sleft/60 << "m" 
						<< (sleft % 60) << "s " 
						<< std::flush;
		}
		// Imprimir el caso especial 100% completado
		if (iEvent == nEvents-1){
			end = std::chrono::steady_clock::now();
			cout << endl;
			cout << "\r" << fixed  <<setprecision(1) << "100.00%" << std::flush << endl;
		}


	} // Main Loop Ends

	// Guardar tree en ROOT file
	f1.cd();
	Data.Write("",2);

	// Imprimir info relevante (parámetros y tiempo)
	printf("\nMAIN PARAMETERS:\n");
	printf("  nEvents = %sM = %d\n", argv[1], nEvents);
	printf("  seed = %d\n",seed);

	int stotal = std::chrono::duration_cast<std::chrono::seconds>
				(end - start).count();
	cout 	<< endl << "Tiempo transcurrido: " 
				<< stotal/3600 << "h" 
				<< stotal/60 << "m" 
				<< (stotal % 60) << "s" << endl;

	// End of file
	cout << endl << "SUCCESS!" << endl;


}