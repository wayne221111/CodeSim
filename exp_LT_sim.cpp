/*
 *  exp_LT_sim.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/11/14.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>    
#include <cstdlib>
#include <vector>
#include <math.h>
#include "randomc.h"
#include "LT.h"
#include <algorithm>
#include <omp.h>


#define K 1000

#define RUN 1000
#define C 0.05
#define Delta 0.02
#define STEPS 16
#define MaxN (K*(1+Delta*(STEPS-1)))

using namespace std;
using namespace CodeSim;

//int M_encode[MaxN][K];				// 記錄 完整的 matrix 
//int V_recover[K];					// 記錄 哪些 symbols 已經被 recover
//int S_recover = K;					// 記錄 有幾個 bit 已解碼  = sum(V_recover);
//int V_degree[MaxN];					// 記錄 M_encode 每個 code word 還有多少 degrees
//int V_ripp[K];						// 記錄 還沒被使用的 degree one symbol    
//int S_ripp = 0;						// 記錄 V_ripp 的個數 = nnz(V_ripp);    
//
//int M_code2sym[MaxN][K];			// 記錄每個 code 哪些 bit 有值
//int S_code2sym[MaxN];
//int M_sym2code[K][MaxN];            // 記錄反向連結 - 每個symbol連到哪幾個 codeword
//int S_sym2code[K];                  // 反向連結的 size 

//unsigned long ErrorCount[STEPS][16];
double BER[STEPS][RUN];
int Dsize = 10;

int		Degree[10] = //{1,2,3,4,5,8,9,19,65,66};
{1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
double  Omega[10] = { 0.060728,0.493283,0.178695,0.033649,0.090168,
	0.008898,0.011950,0.013815, 0.005120,0.103694};

//uniformRandom *Rnd;
int g_seed = (int)time(0);
CRandomMersenne Rnd(g_seed);


double* D;
double* SD;


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


int main(){
	//int i;
	//Rnd = new ran3();
	
	//D = Robust_Soliton_Distribution(K,C,Delta);
	
	//cout << "haha";
	
	D = &Omega[0];
	
	SD = new double[Dsize];
	
	//while (cin >> D[0])
	{
		for (int i =0; i<10; i++) {
			cin >> D[i];
		}
//		for (int i = 0; i<STEPS; i++) {
//			//BER[i] = 0;
//			for(int j = 0; j< 16; j++)
//				ErrorCount[i][j] = 0;
//		}
		
		
		//int n=0;
		#pragma omp parallel for num_threads(6)
		for(int run=0;run<RUN;run++){
			LT_sim<Bit> sim(K, MaxN, Dsize, Degree, D, Rnd.BRandom());
			//cout << "RUN " << i <<'\n';
			//#pragma omp parallel for
			for (int i = 0; i< STEPS; i++) {
				sim.seqReceive( K*(1+Delta*i) -1);
				//sim.decode();
				double t = sim.failureRate();// = Encoder(K, K*(1.05+0.01*i), Dsize);
//				#pragma omp atomic
//				ErrorCount[i][(int)(t*16)]++;
//				#pragma omp atomic
				BER[i][run]=t;
			}
			
			//cout << '\n';
		}
//		cout<<"Histogram:\n";
//		for (int i = 0; i<16; i++) {
//			for(int j = 0; j< STEPS; j++)
//				cout <<  ErrorCount[j][i] / (double)RUN<< ' ';
//			cout << '\n';
//		}
//		
//		cout << "BER:\n";
//		#pragma omp parallel for num_threads(6)
//		for (int i = 0; i<STEPS; i++) {
//			sort(BER[i], BER[i]+RUN);
//			//cout <<  BER[i][]<< ' ';
//			
//		}
		
		for (int i=0; i<RUN; i++) {
			for (int j=0; j<STEPS; j++) {
				cout << BER[j][i] << '\t';
			}
			cout << '\n';
		}
			
	}
	
	return 0;
}




