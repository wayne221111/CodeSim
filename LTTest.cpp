/*
 *  LTTest.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/2.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "LTCode.h"
#include "randomc.h"
#include<ctime>
#define L 300000
#define RUN 30
using namespace CodeSim;

int main(){
	int Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
	double Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
		5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
	double sum[16];
	for (int i = 0; i< 16; i++)
		sum[i] = 0;
	
	CRandomMersenne r(time(0));
	int startTime = time(0);
	int K = 1000;
	for (int run=0; run<RUN; run++) {
		//cout << "Run: " << run << "\n";
		LTCode<Bit> lt(K, K*1.2, 10, Degree_of_Edge, Omega, r.BRandom());
		
		Codeword<Bit> in(L, 0);
		for (int i=0; i<in.size(); i+=2) {
			in[i] = r.IRandomX(0, 1);
		}
		Codeword<Bit> c = lt.encode(in);
		
		for (int q=0; q*K<L; q++){
			for (int i=K*1.04; i<K*1.2; i++) {
				c[q*K*1.2 + i].setErased(true);
			}
		}
		
		for (int i = 0; i< 16; i++) {
			for (int q=0; q*K<L; q++){
				for (int j=K*(1.04+0.01*i); j<K*(1.05+0.01*i); j++) {
					c[q*K*1.2 +j].setErased(false);
				}
			}
			Codeword<Bit> out = lt.decode(c);
			
			
			//sum[i] += sim.run();
			//cout << sim.run();
			int s=0;
			for (int j = 0; j<out.size(); j++) {
				if ( !out[j].isErased() && in[j] == out[j]) 
				{
					s++;
				}
				
			}
			sum[i] += 1- s/(double)( L) ;
		}
	}
	for (int i = 0; i< 16; i++)
		cout << sum[i] / RUN  << '\t';
	
	cout << "Time: " << time(0) - startTime << '\n';
	
    return 0;
	
}