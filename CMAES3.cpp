/*
 *  main_burst_steko.cpp
 *  LT_CMAES
 *
 *  Created by 刁培倫 on 2010/7/21.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "cmaes_interface.h"
#include "randomc.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>    
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <climits>
#include <cmath>
#include <string>

#include <omp.h>
#include "LT.h"
#include "DynamicFitting.h"

// =================================================================
//#define K 1000			// K size
//#define MaxN 10500		// set Max code word number (set 2*K ,but we only use 1.2*K)
//#define Run 10000			// how many simulations per fitness 
#define MAXFEC 10000		// set Max function evaluations in CMAES  
#define Lambda 10		// set parameter lambda in CMAES
#define INFO 1			// 1 : show the info during evolution , 0 : don't display
// =================================================================
//#define Delta 0.005
//#define STEPS 61
#define MaxEpsilon (Delta*(STEPS-1))
#define MaxN (K*(1+MaxEpsilon))
#define P_e_min (1/(double)(K*Run))   // minimum of BER
using namespace std;
using namespace CodeSim;

int g_seed = (int)time(0);
CRandomMersenne RanGen(g_seed);
int K, cubic_fitting, error_exponent, STEPS, Run;
double Delta;
// =================================================================
// 依照要跑的 degree 去設定 
int 	Dsize = 10;
//int		Set_tags[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
// =================================================================

int*	Tags;
double *D;
//[] = {5.5866E-03, 5.0000E-01, 1.6667E-01, 8.3333E-02, 5.0000E-02,
//	2.3810E-02, 1.3889E-02, 2.6316E-03, 2.9200E-04, 1.5379E-01};
double* SD;
double* Std;


//#define epsilonIndex 2
//double epsilons[epsilonIndex] = {0.05, 0.15},//{0.05, 0.06, 0.07, 0.08, 0.09, 

//	0.1, 0.11, 0.12, 0.13, 0.14, 
//	0.15, 0.16, 0.17, 0.18, 0.19,
//	0.2},//{0.05, 0.06,0.08, 0.12,0.20},
//targetErrorRate[epsilonIndex] = {0.1, 0.01};//{0.376629,
vector<double> epsilons, epsilons_w, targetErrorRate, targetErrorRate_w;
//	0.349958, 0.237445, 0.18919, 0.108815, 0.0701708,
//	0.0347487, 0.01925, 0.0085005, 0.0044277, 0.0020448,
//	0.0010388, 0.0004038, 0.0003112, 0.0002659, 0.0001987};
//errorRateBound[epsilonIndex]={4,4,4,4,4},
//epsilonBurstBound[epsilonIndex] = {0.5,0.4,0.3,0.2,0.1},
double epsilonBurstBound = 0.5,
errorDensityBound = 0.33333;
//errorDensityBound[epsilonIndex] = {0.33333, 0.26666, 0.2, 0.166666, 0.133333},
//areaWeight[epsilonIndex] = {0, 500, 1000, 2000, 4000};
int winSize = 30;
//int winSize[epsilonIndex] = {30, 30, 30, 30, 30};
//double failurePenalty[epsilonIndex] = {20, 40, 60, 80, 100}
;







istream& mygetline ( istream& is, string& str );

double* normolize(double* d){
	int i;
	double z = 0;
	for(i=0;i<Dsize;i++){
		if(d[i]<0) d[i] = -d[i];
		z = z + d[i];
	}
	for(i=0;i<Dsize;i++) d[i] = d[i]/z;
	return d;
}

//  初始化設定參數  從uniform distribution 開始  STD 設為 0.025 
void Parameter_init(){
	//Tags = new int[Dsize];
	//D  = new double[Dsize];
	
	SD  = new double[Dsize];
	Std    = new double[Dsize];
	for(int i =0;i<Dsize;i++){
		//Tags[i] = Set_tags[i];
		//D[i] = 1/(double)Dsize;
		Std[i] 	 = 0.1;
	}
}

double weighting(double epsilon, double f_rate){
	return epsilon*pow(f_rate,1.5);
	
}


inline double exceed_penalty(double value, double base, double penalty_ratio)
{
	if (value <= base) {
		return 0;
	}
	
	return (value - base) * penalty_ratio;
}

double fitfun(double* Indiv , int dim, bool &needResample, vector<double> &parameters){
	
	normolize(Indiv);
	
	parameters.assign(epsilons.size()+targetErrorRate.size(),0);
	
	double fit=0, *err = new double[STEPS];
	unsigned long	*errorCount = new unsigned long[STEPS];
	double *failureCount = new double[STEPS];
	for (int i=0; i<STEPS; i++) {
		err[i]=0;
		failureCount[i] = 0;
		errorCount[i] = 0;
	}
	
	// run simulation
	#pragma omp parallel for schedule(dynamic) num_threads(6) reduction(+:fit)
	for(int i=0;i<Run;i++){
		
		LT_sim<Bit> sim(K, (int) MaxN, Dsize, Tags, Indiv, RanGen.BRandom());
		
		for (int i=0; i<STEPS; i++) {
			sim.seqReceive(K*(1+Delta*(i))-1);
			//sim.decode();
			double temp = sim.failureRate();
			#pragma omp atomic
			errorCount[i] += K*temp;
			
			
		}
		
	}//for
	
	vector<double> x(STEPS, 0.0);
	
	for (int i=0; i<STEPS; i++) {
//		if (failureCount[i] > 80) {
//			//fit +=	200;
//			//needResample = true;
//		}
		if (errorCount[i] > 0) {
			err[i] = errorCount[i]/(double)(K*Run);
			err[i] = log10(err[i]);
			
			if(error_exponent)
				err[i] /= 1+i*Delta;
		}
		else {
			err[i] = log10(P_e_min);
		}
		
		// prepare x
		x[i] = i*Delta;
		
	}
	
	if(cubic_fitting){
		// start regression
		vector<double> a(4, 0.0), y(err, err+STEPS);
		double ss;
		
		regression(x, y, a, cubic, ss);
		

		


			
		for (int i=0; i<epsilons.size(); i++) {
			fit += (eval(epsilons[i], a, cubic)-log10(P_e_min)) / abs(log10(P_e_min)) * epsilons_w[i];
			parameters[i] = (eval(epsilons[i], a, cubic)-log10(P_e_min)) / abs(log10(P_e_min));
		}
			
		for (int i=0; i<targetErrorRate.size(); i++) {	
			
			vector<double> a2 = a;
			a2[0] -= log10(targetErrorRate[i]);
			vector<double> sols = find_sol(a2, 0, (STEPS-1)*Delta, STEPS);
			double solution;
			switch (sols.size()) {
				case 0:
					if (eval(0, a, cubic) < log10(targetErrorRate[i])) {
						solution = 0;
					}
					else {
						solution = (STEPS-1)*Delta;
					}
					break;
				case 1:
					solution = sols[0];
					break;
				

				default:
					solution = *sols.rbegin();
					break;
			}
			
			fit += solution/(Delta*(STEPS-1)) * targetErrorRate_w[i];
			parameters[epsilons.size()+i] = solution/(Delta*(STEPS-1));
		}
		
	}
	
	else {
		for (int i=0; i<epsilons.size(); i++) {
			fit += (err[(int)(epsilons[i]/Delta)]-log10(P_e_min))  / abs(log10(P_e_min)) * epsilons_w[i];
			parameters[i] = (err[(int)(epsilons[i]/Delta)]-log10(P_e_min))  / abs(log10(P_e_min));
		}
		
		for (int i=0; i<targetErrorRate.size(); i++) {	
			
			double min_diff=999;
			int min_i;
			for (int j=0; j<STEPS; j++) {
				if ( abs(err[j]-log10(targetErrorRate[i])) < min_diff) {
					min_diff = abs(err[j]-log10(targetErrorRate[i]));
					min_i = j;
				}
			}
			

			fit += min_i/(double)(STEPS-1) * targetErrorRate_w[i];
			parameters[epsilons.size()+i] = min_i/(double)(STEPS-1);
		}
	}

	
	delete [] err;
	delete [] errorCount;
	delete [] failureCount;
	return fit;
}

/* the optimization loop */
int main(int argn, char **args) {	
	ifstream ifs;
	string filename;
	if (argn == 1) {
		cerr << "Usage: CMAES3.out <filename>" << endl;
		//exit(1);
		filename = "input_script";
		ifs.open(filename.c_str());
	}
	else {
		ifs.open(args[1]);
		filename = args[1];
		if(ifs.fail()){
			cerr << "Error: can not open file \"" << args[1] << '\"' << endl;
			exit(1);
		}
	}
	string comm;
	
	// read K
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		iss >> K >> Run;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	// read Number of tags
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		iss >> Dsize;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}

	
	// read Tags
	Tags = new int[Dsize];
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		for (int i=0; i<Dsize; i++) {
			iss >> Tags[i];
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	
	// read initial distribution
	D  = new double[Dsize];
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		for (int i=0; i<Dsize; i++) {
			iss >> D[i];
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read STEPS and Delta
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		iss >> STEPS >> Delta;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read cubic fitting
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		iss >> cubic_fitting;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}	
	
	// read error exponent
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		iss >> error_exponent;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}	
	
	// read epsilons
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		double t;
		while (iss >> t) {
			epsilons.push_back(t);
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read epsilons weighting
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		double t;
		while (iss >> t) {
			epsilons_w.push_back(t);
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	if(epsilons.size() != epsilons_w.size()){
		cerr << "Error: epsilons list size doesn't match to weighting list" << endl;
		exit(1);
	}
	
	// read targetErrorRate
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		double t;
		while (iss >> t) {
			targetErrorRate.push_back(t);
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read targetErrorRate weighting
	if(mygetline(ifs,comm)){
		istringstream iss(comm);
		double t;
		while (iss >> t) {
			targetErrorRate_w.push_back(t);
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	if(targetErrorRate.size() != targetErrorRate_w.size()){
		cerr << "Error: targetErrorRate list size doesn't match to weighting list" << endl;
		exit(1);
	}
	
	
	int i; 
	fstream fs;	
	//Rnd = new ran0(0);
  	time_t rawtime;
  	struct tm * timeinfo;
	
	cmaes_t* evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xbest;
	
	// open file
	//cout << "Enter filename: ";
	comm = "result_";
	comm += filename;
	//getline(cin, comm);
	fs.open(comm.c_str(),fstream::out);
	
	// recored start time 
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"Start time : "<<asctime(timeinfo)<<endl;	
	
//	cout << "Enter Comment: ";
//	getline(cin, comm);
	
	Parameter_init();
	// write Tags and init distribtuion into file
	fs << "Comment: " << comm << '\n';
	fs<<"Tags\n";
	for(i=0;i<Dsize;i++) fs<<Tags[i]<<"\t";
	fs<<"\nInitial distribution \n";
	for(i=0;i<Dsize;i++) fs<<D[i]<<"\t";
	fs<<"\nEpsilons \n";
	for(i=0;i<epsilons.size();i++) fs<<epsilons[i]<<"\t";
	fs<<"\nEpsilons weighting\n";
	for(i=0;i<epsilons_w.size();i++) fs<<epsilons_w[i]<<"\t";
	fs<<"\nTarget Error Rate \n";
	for(i=0;i<targetErrorRate.size();i++) fs<<targetErrorRate[i]<<"\t";
	fs<<"\nTarget Error Rate weighting\n";
	for(i=0;i<targetErrorRate_w.size();i++) fs<<targetErrorRate_w[i]<<"\t";
	fs<<"\nGen\tFEvals\tFitness\tFbest\tXbest dist.\t";
	for(i=0;i<Dsize;i++) fs<< Tags[i] << '\t';
	fs<< "para.\t";
	for(i=0;i<epsilons.size();i++) fs<<"p@e="<<epsilons[i]<<"\t";
	for(i=0;i<targetErrorRate.size();i++) fs<<"e@p="<<targetErrorRate[i]<<"\t";
	fs<<endl;
	
	evo = new cmaes_t();
	
	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(evo, Dsize, D, Std, 0, Lambda, "non"); 
	evo->sp.stopMaxFunEvals = MAXFEC;	
	cout<<cmaes_SayHello(evo)<<endl;
	
	/* Iterate until stop criterion holds */
	while(!cmaes_TestForTermination(evo))
	{ 						
		vector<double> xbest_parameters;
		double min_fit = 9999;
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(evo); /* do not change content of pop */
		/* evaluate the new search points using fitfun from above */ 
		for (i = 0; i < Lambda; ++i) {
			bool needResample = false;
			vector<double> t_parameters;
			arFunvals[i] = fitfun(pop[i], Dsize, needResample, t_parameters);
			
			if (needResample) {
				pop = cmaes_ReSampleSingle(evo, i);
				i--;
				cout << "R";
			}
			else {
				cout << "E";
			}
			cout.flush();
			if(min_fit > arFunvals[i]){
				min_fit = arFunvals[i];
				xbest_parameters = t_parameters;
			}
		}
		cout << '\n';
		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(evo, arFunvals);  
		/* read instructions for printing output or changing termination conditions */ 
		xbest = normolize(cmaes_GetNew(evo, "xbest"));
		fs<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<"\t";
		fs.setf(ios::fixed);
		fs.precision(6);
		fs << "dist.\t";
		for(i=0;i<Dsize;i++) 
			fs<<setw(8)<<xbest[i]<<"\t";
		
		fs << "para.\t";
		for (int i=0; i<xbest_parameters.size(); i++) {
			fs<<setw(8)<<xbest_parameters[i]<<"\t";
		}
		
		fs.unsetf(ios::fixed);
		fs<<endl;
		
		if(INFO==1) cout<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<endl;
		//fflush(stdout); /* useful in MinGW */
		
		string dist = "dist_"+filename;
		ofstream dd(dist.c_str());
		dd << K << endl;
		dd << Dsize << endl;
		for(int i=0;i<Dsize;i++) 
			dd<<Tags[i]<<"\t";
		dd << endl;
		for(int i=0;i<Dsize;i++) 
			dd<<setw(12)<<xbest[i]<<"\t";
		dd << endl;
		if(K==10000)
			dd << "16 0.013" << endl;
		else {
			dd << "16 0.03" << endl;
		}

		dd.close();
		delete [] xbest;
	}
	
	printf("Stop:\n%s\n",  cmaes_TestForTermination(evo)); /* print termination reason */
	cmaes_WriteToFile(evo, "all", "allcmaes.dat");         /* write final results */
	
	// recored end time 		
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"\nStop time : "<<asctime(timeinfo)<<endl;
	fs.close();	
	
	cmaes_exit(evo); /* release memory */ 		
	delete [] Tags;
	
	delete [] SD;
	delete [] Std;	
	delete (evo);
	
	
	cout << "Running histogram...";
	cout.flush();
	string cmd = "~/CodeSim/histogram.out < dist_" + filename 
										+ " > histo_" + filename;
	
	system(cmd.c_str() );
	cout << "done!" << endl;
	return 0;
}

istream& mygetline ( istream& is, string& str ){
	string s;
	while(getline(is, s)){
		if (s.size() == 0 || s[0] != '#') {
			str = s;
			break;
		}
	}
	
	return is;
}
