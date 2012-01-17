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

// =================================================================
#define K 1000			// K size
//#define MaxN 10500		// set Max code word number (set 2*K ,but we only use 1.2*K)
#define Run 10000			// how many simulations per fitness 
#define MAXFEC 10000		// set Max function evaluations in CMAES  
#define Lambda 20		// set parameter lambda in CMAES
#define INFO 1			// 1 : show the info during evolution , 0 : don't display
// =================================================================
#define Delta 0.001
#define STEPS 501
#define MaxN (K*(1+Delta*(STEPS-1)))
using namespace std;
using namespace CodeSim;

int g_seed = (int)time(0);
CRandomMersenne RanGen(g_seed);

// =================================================================
// 依照要跑的 degree 去設定 
int 	Dsize = 10;
int		Set_tags[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
// =================================================================

int*	Tags;
double D[] = {5.5866E-03, 5.0000E-01, 1.6667E-01, 8.3333E-02, 5.0000E-02,
	2.3810E-02, 1.3889E-02, 2.6316E-03, 2.9200E-04, 1.5379E-01};
double* SD;
double* Std;


#define epsilonIndex 2
double epsilons[epsilonIndex] = {0.05, 0.15},//{0.05, 0.06, 0.07, 0.08, 0.09, 
//	0.1, 0.11, 0.12, 0.13, 0.14, 
//	0.15, 0.16, 0.17, 0.18, 0.19,
//	0.2},//{0.05, 0.06,0.08, 0.12,0.20},
targetErrorRate[epsilonIndex] = {0.1, 0.001};//{0.376629,
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






double e = 1.05;					

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
	Tags = new int[Dsize];
	//D  = new double[Dsize];
	
	SD  = new double[Dsize];
	Std    = new double[Dsize];
	for(int i =0;i<Dsize;i++){
		Tags[i] = Set_tags[i];
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

double fitfun(double* Indiv , int dim, bool &needResample){
	
	normolize(Indiv);
	
	double fit=0, err[STEPS];
	long	errorCount[STEPS];
	double failureCount[STEPS];
	for (int i=0; i<STEPS; i++) {
		err[i]=0;
		failureCount[i] = 0;
		errorCount[i] = 0;
	}
	
	
#pragma omp parallel for num_threads(6) reduction(+:fit)
	for(int i=0;i<Run;i++){
		//cout << "Run "<< i+1 << endl;
		//Codeword<Bit> decodePattern[epsilonIndex];
		
		{
			LT_sim<Bit> sim(K, (int) MaxN, Dsize, Set_tags, Indiv, RanGen.BRandom());
			
			for (int i=0; i<STEPS; i++) {
				sim.seqReceive(K*(1+Delta*(i))-1);
				//sim.decode();
				double temp = sim.failureRate();
				#pragma omp atomic
				errorCount[i] += K*temp;
				
				//				if (temp > epsilonBurstBound) {
				//					#pragma omp atomic
				//					failureCount[i] += 1;
				//				}
				//				Codeword<Bit> t = sim.getResult();
				//				
				//				decodePattern[i].insert(decodePattern[i].end(), t.begin(), t.end());
				
			}
		}
		
		
		//		for (int i=0; i<epsilonIndex; i++) {
		//			
		//			int errNO=0, errLen=0;
		//			for (int p=0; p<winSize; p++) {
		//				if (decodePattern[i][p].isErased()) {
		//					errNO ++;
		//				}
		//			}
		//			if(errNO/(double)winSize > errorDensityBound)
		//				errLen=1;
		//			
		//			for (int p=winSize; p< 100*K; p++) {
		//				if (decodePattern[i][p].isErased()) {
		//					errNO++;
		//				}
		//				if (decodePattern[i][p-winSize].isErased()) {
		//					errNO --;
		//				}
		//				
		//				if(errNO/(double)winSize > errorDensityBound)
		//				{
		//					errLen ++;
		//				}
		//				else {
		//					if (errLen > 750) {
		//						//fit +=failurePenalty[i];
		//						#pragma omp atomic
		//						failureCount[i] += errLen / 750.0;
		//					}
		//					errLen = 0;
		//				}
		//				
		//				
		//			}	
		//			
		//			
		//		}
		
	}
	
	
	
	for (int i=0; i<STEPS; i++) {
		if (failureCount[i] > 80) {
			//fit +=	200;
			needResample = true;
		}
		if (errorCount[i] > 0) {
			err[i] = errorCount[i]/(double)(K*Run);
			err[i] = log10(err[i]);//+4;
			//err[i] += exceed_penalty(err[i], errorRateBound[i], 1000);
			//fit += abs ( err[i] - log10(targetErrorRate[i]) );
		}
		else {
			err[i] = -4;
		}

	}
	
	for (int i=0; i<epsilonIndex; i++) {
		if(err[(int)(epsilons[i] / Delta)] > log10(targetErrorRate[i]))
			needResample = true;
		
		fit += err[(int)(epsilons[i] / Delta)]+4;
		
		double min_diff=999;
		int min_i;
		for (int j=0; j<STEPS; j++) {
			if ( abs(err[j]-log10(targetErrorRate[i])) < min_diff) {
				min_diff = abs(err[j]-log10(targetErrorRate[i]));
				min_i = j;
			}
		}
		
		if (min_i*Delta > epsilons[i]) {
			needResample = true;
		}
		fit += min_i/100.0;
	}
	
	
	
	//delete [] err;
	return fit;
}

/* the optimization loop */
int main(int argn, char **args) {	
	
	int i; 
	fstream fs;	
	//Rnd = new ran0(0);
  	time_t rawtime;
  	struct tm * timeinfo;
	
	cmaes_t* evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xbest;
	
	// open file
	cout << "Enter filename: ";
	string comm;
	getline(cin, comm);
	fs.open(comm.c_str(),fstream::out);
	
	// recored start time 
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"Start time : "<<asctime(timeinfo)<<endl;	
	
	cout << "Enter Comment: ";
	getline(cin, comm);
	
	Parameter_init();
	// write Tags and init distribtuion into file
	fs << "Comment: " << comm << '\n';
	fs<<"Tags\n";
	for(i=0;i<Dsize;i++) fs<<Tags[i]<<"\t";
	fs<<"\nInitial distribution \n";
	for(i=0;i<Dsize;i++) fs<<D[i]<<"\t";
	fs<<"\nEpsilons \n";
	for(i=0;i<epsilonIndex;i++) fs<<epsilons[i]<<"\t";
	fs<<"\nTarget Error Rate \n";
	for(i=0;i<epsilonIndex;i++) fs<<targetErrorRate[i]<<"\t";
	fs<<"\nGen\tFEvals\tFitness\tFbest\tXbest\n";
	
	evo = new cmaes_t();
	
	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(evo, Dsize, D, Std, 0, Lambda, "non"); 
	evo->sp.stopMaxFunEvals = MAXFEC;	
	cout<<cmaes_SayHello(evo)<<endl;
	
	/* Iterate until stop criterion holds */
	while(!cmaes_TestForTermination(evo))
	{ 						
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(evo); /* do not change content of pop */
		/* evaluate the new search points using fitfun from above */ 
		for (i = 0; i < Lambda; ++i) {
			bool needResample = false;
			arFunvals[i] = fitfun(pop[i], Dsize, needResample);
			if (needResample) {
				pop = cmaes_ReSampleSingle(evo, i);
				i--;
				cout << "R";
			}
			else {
				cout << "E";
			}
			cout.flush();
			
		}
		cout << '\n';
		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(evo, arFunvals);  
		/* read instructions for printing output or changing termination conditions */ 
		xbest = normolize(cmaes_GetNew(evo, "xbest"));
		fs<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<"\t";
		fs.setf(ios::fixed);
		fs.precision(6);
		for(i=0;i<Dsize;i++) 
			fs<<setw(8)<<xbest[i]<<"\t";
		fs.unsetf(ios::fixed);
		fs<<endl;
		
		if(INFO==1) cout<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<endl;
		//fflush(stdout); /* useful in MinGW */
		
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
	delete Tags;
	//delete D;
	delete SD;
	delete Std;	
	delete(evo);
	//delete(Rnd);
	system("pause");
	return 0;
}


