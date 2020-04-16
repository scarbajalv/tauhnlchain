// ./dsprodcppv2 ebeam index

#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <iomanip>
#include <math.h>
#include <string>
#include <cmath>

using namespace std;

double cdfsigmax(double x, double n){
	double cdfsigmaxreturn;
	if(x>=0 && x<0.5){
		cdfsigmaxreturn = 1.*pow(2*x,1/(n+1))-1;
	}
	else{
		cdfsigmaxreturn = 1-1.*pow(2-2*x,1/(n+1));
	}
	return cdfsigmaxreturn;
}
double cdfsigmapt2(double x, double b){
	double sigmapt2return = (1./b)*log(1./(1-x));
	return sigmapt2return;
}


int main(int argc, char *argv[]){

	if (argc <  2 ){
		cout << "ERROR: " <<  endl << "./dsprodcppvX index" << endl;
		exit(EXIT_FAILURE);
	}
	
	// AQUÍ LOS PRINCIPALES PARÁMETROS
	double ebeam = 120.;
	int nevents = 1000;
	double m = 1.96834; // mass of Dmeson
	double a = 15.2431907; // x-gauss
	double b = 1.08; // pT2-exp
		
	//double ebeam = stod(argv[1]);
	double sigma_param = 1./sqrt(2*a);
  default_random_engine generator;
  normal_distribution <double> distribution(0.0, sigma_param);
  double number;

	vector <double> row;
	vector < vector <double> > ivector;
	vector < vector <double> > ovector;

	srand (time(NULL));

	double mbeam = 0.968;
	double mtarget = 0.968;

	double sqrts = sqrt(pow(ebeam+mtarget,2)-(pow(ebeam,2)-pow(mbeam,2)));
	double pzcm, pt, ecm, gamma, v, pz, e, theta, phi, px, py;
	double x_aux pT2_aux;



	// Generamos la rawdata
	cout<<"Generando rawdata..."<<endl;
	for (int i=0; i<nevents; ++i){
		row.clear();
		row.push_back(cdfsigmax( (rand() % 10000001)/10000000. , 6.1));
		number = distribution(generator);
		row.push_back(number);
		ivector.push_back(row);
		cout<<"\r"<<i+1;
	}
	
	cout<<endl;

	// TRANSFORMAR RAWDATA A DSDATA

	
	cout<<"Generando dsdata..."<<endl;
	for (int i=0; i<ivector.size(); ++i){
		row.clear();
		pzcm = 0.5*sqrts*ivector[i][0];
		pt = sqrt(1.*ivector[i][1]);
		ecm = sqrt(1.*pow(pzcm,2)+pow(pt,2)+pow(m,2));
		gamma = (0.5*sqrts)/mtarget;
		v = sqrt(1-1./pow(gamma,2));
		pz = 1.*gamma*(pzcm-(-v)*ecm);
		e = sqrt(1.*pow(pz,2)+pow(pt,2)+pow(m,2));
		theta = atan(pt/pz);
		phi = (rand() % 1001)/1000.*2*3.14159265358;
		//px = sqrt(pow(pz,2)+pow(pt,2))*sin(theta)*cos(phi);
		//py = sqrt(pow(pz,2)+pow(pt,2))*sin(theta)*sin(phi);
		// Debug en caso de que cdfsigmapt2 se haga infinito (x=1)
		if (isnan(e)||isnan(theta)||isinf(e)||isinf(theta)){
			cout<<"nan or inf found!"<<endl;
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

	stringstream ofilenamess;
	ofilenamess<<"dsdata-"<<ebeam<<"-"<<argv[2]<<".dat";
	string ofilename=ofilenamess.str();
	ofstream ofile;
	ofile.open(ofilename);

	cout<<"Exportando dsdata..."<<endl;
	for (int i=0; i<ovector.size(); ++i){
		for (int j=0; j<3; ++j){
			ofile<<setprecision(10)<<ovector[i][j]<<" ";
		}
		ofile<<endl;
		cout<<"\r"<<i+1;
	}

	cout<<endl<<"SUCCESS!"<<endl;
	
	return 0;
}