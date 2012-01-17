#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>    
#include <cstdlib>
#include <vector>
#include <math.h>
#include "randomc.h"
#include "LT.h"
#include "statistics.h"
#include <omp.h>
#include <algorithm>

int K;

int Run;
//#define C 0.05
double Delta;
int STEPS;
#define MaxN (K*(1+Delta*(STEPS-1)))
double targetRho, targetFailureRate, targetEpsilon;
int optimParameter;

using namespace std;
using namespace CodeSim;


//unsigned long ErrorCount[STEPS][16];
//double BER[STEPS];
//vector< unsigned long > BER;
int Dsize;


int* Tags;
//double* Omega;
//int		Degree[10] = //{1,2,3,4,5,8,9,19,65,66};
//		{1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
//double  Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02,
//				5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02
//};

//uniformRandom *Rnd;
int g_seed = (int)time(0);
CRandomMersenne RanGen(g_seed);


double* Indiv;
double* SD;


double e = 1.05;					



int main(){

	cin >> K;
	cin >> Run;
	cin >> Dsize;
	Tags = new int[Dsize];
	Indiv = new double[Dsize];
	
	for (int i=0; i<Dsize; i++) {
		cin >> Tags[i];
	}

	cin >> STEPS >> Delta >> targetRho >> targetFailureRate >> targetEpsilon >> optimParameter;
	
	
	
	while (cin >> Indiv[0])
	{
		for (int i =1; i<Dsize; i++) {
			cin >> Indiv[i];
		}
		
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
//					if(t == 1)
//						needResample = true;
					cout << log10(t);
				}
				else {
					cout << log10(1.0/K) - 0.1;
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
//					if(t == 1)
//						needResample = true;
					cout << log10(t);
				}
				else {
					cout << log10(1.0/Run) -0.1;
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
				
				
				double output = -1;
				for (int i=0; i<STEPS; i++) {
					if (errorCount[i] <= Run*targetFailureRate) {
						output = i*Delta;
						break;
					}
				}
				
				// can not find acceptable epsilon, maximum returned.
				if (output == -1) {
					output =  STEPS*Delta;
				}
				
				cout << output;
			}
				break;
		}
		
		
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  BER[i] / (double)(K*Run)<< '\t';
//		}
		cout << endl;
		
	}
	
	return 0;
}




