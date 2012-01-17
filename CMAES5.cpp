/*
 *  CMAES5.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/10/13.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
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
#include <sstream>
#include <algorithm>
#include <unistd.h>
// =================================================================
//#define K 1000			// K size
//#define MaxN 10500		// set Max code word number (set 2*K ,but we only use 1.2*K)
//#define Run 10000			// how many simulations per fitness 
#define MAXFEC 10000		// set Max function evaluations in CMAES  
#define Lambda 20		// set parameter lambda in CMAES
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
int K, STEPS, Run;
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



double targetRho, targetFailureRate, targetEpsilon;
int optimParameter;

double epsilonBurstBound = 0.5,
errorDensityBound = 0.33333;
int winSize = 30;


istream& mygetline ( istream& is, string& str );

//double* normolize(double* d){
//	int i;
//	double z = 0;
//	for(i=0;i<Dsize;i++){
//		if(d[i]<0) d[i] = -d[i];
//		z = z + d[i];
//	}
//	for(i=0;i<Dsize;i++) d[i] = d[i]/z;
//	return d;
//}

//  初始化設定參數  從uniform distribution 開始  STD 設為 0.025 
void Parameter_init(){
	//Tags = new int[Dsize];
	//D  = new double[Dsize];
	
	SD  = new double[Dsize];
	Std    = new double[Dsize-1+Dsize-3];
	for(int i =0;i<Dsize-1+Dsize-3;i++){
		//Tags[i] = Set_tags[i];
		//D[i] = 1/(double)Dsize;
		if(i<Dsize-1) {
			Std[i] 	 = 0.025;
		}
		else {
			Std[i] = 10;
		}

	}
	
}


inline double exceed_penalty(double value, double base, double penalty_ratio)
{
	if (value <= base) {
		return 0;
	}
	
	return (value - base) * penalty_ratio;
}

double fitfun(double* Indiv , int dim, bool &needResample){
	
	
	switch (optimParameter) {
		case 1:// optimize rho
		{
			vector<double> failureRatios(Run, 0);
			// run simulation
			#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
			for(int run=0;run<Run;run++){
				
				int seed;
				#pragma omp critical
				seed = RanGen.BRandom();
				
				LT_sim<Bit> sim(K, K*(1+targetEpsilon)+0.1, Dsize, Tags, Indiv, seed);
				
				
				sim.seqReceive(K*(1+targetEpsilon)-0.9);
				failureRatios[run] = sim.failureRate();
				
				
			}// end of simulation
			
			sort(failureRatios.begin(), failureRatios.end());
			
			double t = failureRatios[Run*(1-targetFailureRate)-0.9];
			if(t > 0) {
				if(t == 1)
					needResample = true;
				return log10(t);
			}
			else {
				return log10(1.0/K) - 0.1;
			}
			
			
		}
			break;
		case 2:// optimize p
		{
			vector<double> failureRatios(Run, 0);
			// run simulation
			#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
			for(int run=0;run<Run;run++){
				
				int seed;
				#pragma omp critical
				seed = RanGen.BRandom();
				
				LT_sim<Bit> sim(K, K*(1+targetEpsilon)+0.1, Dsize, Tags, Indiv, seed);
				
				
				sim.seqReceive(K*(1+targetEpsilon)-0.9);
				failureRatios[run] = sim.failureRate();
				
				
			}// end of simulation
			
			sort(failureRatios.begin(), failureRatios.end());
			
			double t = 1;
			for (int i=failureRatios.size()-1; i>=0; i--) {
				if(failureRatios[i] <= targetRho) {
					t= (Run-i-1)/(double)Run;
					break;
				}
			}
			
			if(t > 0) {
				if(t == 1)
					needResample = true;
				return log10(t);
			}
			else {
				return log10(1.0/Run) -0.1;
			}
			
			
		}
			break;
		case 3: // optimize epsilon
		{
			vector<int> errorCount(STEPS, 0);
			
			// run simulation
			#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
			for(int run=0;run<Run;run++){
				
				int seed;
				#pragma omp critical
				seed = RanGen.BRandom();
				
				LT_sim<Bit> sim(K, MaxN+0.1, Dsize, Tags, Indiv, seed);
				
				for (int i=0; i<STEPS; i++) {
					sim.seqReceive(K*(1+Delta*(i))-0.9);
					int temp = sim.getNumOfErased();
					if(temp > K*targetRho)
					{
						#pragma omp atomic
						errorCount[i] += 1;
					}
					else {
						break;
					}
				}
				
			}// end of simulation
			
			for (int i=0; i<STEPS; i++) {
				if (errorCount[i] <= Run*targetFailureRate) {
					return i*Delta;
					break;
				}
			}
			
			// can not find acceptable epsilon, maximum returned.
			return STEPS*Delta;
		}
			break;
	}
	
	cerr << "Error!" << endl;
	exit(-1);
	return 0;
	
}

void nearestPoint(double *p) {
	double sum = 0.0, shift = 0.0;
	int positiveDim=0;
	for (int i = 0; i< Dsize-1; i++){
		if(p[i]>0){
			sum += p[i];
			positiveDim = positiveDim + 1;
		}
	}
	if(sum>1) {
		shift = (sum-1) / positiveDim;
	}
	
	for (int i = 0; i< Dsize-1; i++){
		if(p[i]>0){
			p[i] -= shift;
		}
		else {
			p[i] = 0;
		}
	}
	
	if(sum>1)
		nearestPoint(p);
}

/* the optimization loop */
int main(int argn, char **args) {	
	ifstream ifs;
	string filename;
	if (argn == 1) {
		cerr << "Usage: CMAES3.out <filename>" << endl;
		exit(1);
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
	string tmp_string;
	
	// read K
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> K >> Run;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	// read Number of tags
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> Dsize;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	
	// read Tags
	Tags = new int[Dsize];
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
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
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		for (int i=0; i<Dsize; i++) {
			iss >> D[i];
		}
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read STEPS and Delta
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> STEPS >> Delta;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}	
	
	// read minimum acceptable failure ratio rho_tilde
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> targetRho;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read acceptable block failure rate
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> targetFailureRate;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	
	// read targetEpsilon
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> targetEpsilon;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	// read optimization Parameter
	if(mygetline(ifs,tmp_string)){
		istringstream iss(tmp_string);
		iss >> optimParameter;
	}
	else {
		cerr << "inputfile: "<< filename << ": format error"<< endl;
		exit(1);
	}
	
	int i; 
	fstream fs, fit_log;	
	//Rnd = new ran0(0);
  	time_t rawtime;
  	struct tm * timeinfo;
	
	cmaes_t* evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xbest;
	
	
	tmp_string = filename + "_result.txt";
	fs.open(tmp_string.c_str(),fstream::out);
	tmp_string = filename + "_fitness_log.txt";
	fit_log.open(tmp_string.c_str(), fstream::out);
	
	// recored start time 
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"Start time : "<<asctime(timeinfo)<<endl;	
	
	nice(20);
	Parameter_init();
	// write Tags and init distribtuion into file
	fs << "Comment: " << filename << '\n';
	fs<<"Tags\n";
	for(i=0;i<Dsize;i++) fs<<Tags[i]<<"\t";
	fs<<"\nInitial distribution \n";
	for(i=0;i<Dsize;i++) fs<<D[i]<<"\t";
	fs<<"\nminimum acceptable failure ratios \n";
	fs<<targetRho;
	fs<<"\nTarget block failure rate\n";
	fs<<targetFailureRate;
	fs<<"\nTarget Epsilon\n";
	fs<<targetEpsilon;
	fs<<"\nOptimization Parameter\n"<<optimParameter;
	fs<<"\nGen\tFEvals\tFitness\tFbest\tXbest dist.\t";
	for(i=0;i<Dsize;i++) fs<< Tags[i] << '\t';
	fs<<endl;
	
	evo = new cmaes_t();
	
	double DD[Dsize-1+Dsize-3];
	for (int i=0; i<Dsize-1; i++) {
		DD[i] = D[i];
	}
	for (int i=0; i<Dsize-3; i++) {
		DD[i+Dsize-1] = Tags[i+3];
	}
	
	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(evo, Dsize-1+Dsize-3, DD, Std, RanGen.BRandom(), Lambda, "non");
	evo->sp.stopMaxFunEvals = MAXFEC;	
	cout<<cmaes_SayHello(evo)<<endl;
	
	//omp_set_nested(1);
	
	/* Iterate until stop criterion holds */
	while(!cmaes_TestForTermination(evo))
	{ 						
		int resampleTime = 0;
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(evo); /* do not change content of pop */
		/* evaluate the new search points using fitfun from above */ 
		//#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
		for (i = 0; i < Lambda; ++i) {
			bool needResample = true;
			
			
			
			while (needResample) {
				needResample = false;
				
				nearestPoint(pop[i]);
				for (int j=0; j<Dsize-3; j++) {
					if(pop[i][j+Dsize-1]>K) {
						pop[i][j+Dsize-1] = K;
					}
					if(pop[i][j+Dsize-1]<4) {
						pop[i][j+Dsize-1] = 4;
					}
				}
				
				double* dist = new double[Dsize];
				double sum = 0;
				for (int j=0; j<Dsize-1; j++) {
					if (pop[i][j] < 0) {
						needResample = true;
						break;
					}
					dist[j] = pop[i][j];
					sum += pop[i][j];
				}
				dist[Dsize-1] = 1-sum;
				if(dist[Dsize-1] < 0)
					needResample = true;
				
				if(needResample == false) {
					for (int j=0; j<Dsize-3; j++) {
						Tags[j+3] = pop[i][j+Dsize-1]+0.5;
					}
					arFunvals[i] = fitfun(dist, Dsize, needResample);
					if(needResample){
						#pragma omp atomic
						resampleTime++;
						if(resampleTime>1000) {
							needResample = false;
						}
					}
				}
				else {
					arFunvals[i] = 9999;
				}
				
				
				if (needResample) {
					
					pop = cmaes_ReSampleSingle(evo, i);
					
				}
				else {
					cout << "E";
				}
				cout.flush();
				
				delete [] dist;
			}
		}
		cout << '\n';
		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(evo, arFunvals);  
		/* read instructions for printing output or changing termination conditions */ 
		xbest = cmaes_GetNew(evo, "xbest");
		fs<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<"\t";
		fit_log << cmaes_Get(evo, "fitness") << '\t';
		fit_log.flush();
		fs.setf(ios::fixed);
		fs.precision(6);
		fs << "dist.\t";
		double sum = 0;
		for(i=0;i<Dsize-1;i++) { 
			fs<<setw(8)<<xbest[i]<<"\t";
			sum += xbest[i];
		}
		fs<<setw(8)<<1-sum<<"\t";
		for (int j=0; j<Dsize-3; j++) {
			fs<<(xbest[j+Dsize-1]+0.5)<<"\t";
		}
		fs.unsetf(ios::fixed);
		fs<<endl;
		
		if(INFO==1) cout<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<endl;
		//fflush(stdout); /* useful in MinGW */
		
		string dist = filename+"_dist.txt";
		ofstream dd(dist.c_str());
		dd << K << endl;
		dd << Dsize << endl;
		for(int i=0;i<3;i++) 
			dd<<Tags[i]<<"\t";
		for (int j=0; j<Dsize-3; j++) {
			dd<<(int)(xbest[j+Dsize-1]+0.5)<<"\t";
		}
		dd << endl;
		for(int i=0;i<Dsize-1;i++) 
			dd<<setw(12)<<xbest[i]<<"\t";
		dd << setw(12)<< 1-sum << endl;
		//		if(K==10000)
		//			dd << "16 0.013" << endl;
		//		else {
		//			dd << "16 0.03" << endl;
		//		}
		
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
	fit_log.close();
	
	cmaes_exit(evo); /* release memory */ 		
	delete [] Tags;
	
	delete [] SD;
	delete [] Std;	
	delete (evo);
	
	
	cout << "Running histogram...";
	cout.flush();
	string cmd = "./histogram.out " + filename 
	+ "_dist.txt " + filename+ "_histo.txt";
	
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