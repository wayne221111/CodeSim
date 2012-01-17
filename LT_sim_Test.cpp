#include <iostream>
#include <limits>
#include <ctime>
#include "LT.h"
#define Run 20
using namespace std;
using namespace CodeSim;
int main (int argc, char * const argv[]) {
    
	//LT_sim sim;
	int seed = time(0);
	CRandomMersenne r(seed), r2(seed+1);
	
	
	
	int Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
	double Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
				5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
	double sum[16];
	for (int i = 0; i< 16; i++)
		sum[i] = 0;
	
	int K = 1000;
	int startTiem = time(0);
	for (int run=0; run<Run; run++) {
		//cout << "Run: " << run << "\n";
		LT_sim<Bit> sim;
		sim.init(K, K*1.2, 10, Degree_of_Edge, Omega, r.BRandom());
		
		Codeword<Bit> in(K, 0);
		for (int i=1; i<K; i+=2) {
			in[i] = r2.IRandom(0, 2);
		}
		Codeword<Bit> c = sim.encode(in);
		
		for (int i=K*1.04; i<K*1.2; i++) {
			c[i].setErased(true);
		}
		
		for (int i = 0; i< 16; i++) {
//			double t = Encoder(K, K*(1.05+0.01*i), Dsize);
//			ErrorCount[i][(int)(t*16)]++;
//			BER[i]+=t;
			
			
			
//			sim.seqReceive(K*(1.05+0.01*i) );
//			sim.decode();
			
			
			
//			for (int i=0; i< c.size(); i++) {
//				if (r.Random() < 0.1) {
//					c[i].setErased(1);
//				}
//			}
			
			for (int j=K*(1.04+0.01*i); j<K*(1.05+0.01*i); j++) {
				c[j].setErased(false);
			}
			Codeword<Bit> out = sim.decode(c);
			
			
			//sum[i] += sim.run();
			//cout << sim.run();
			int s=0;
			for (int j = 0; j<K; j++) {
				if ( !out[j].isErased() && in[j] == out[j]) 
				{
					s++;
				}
				
			}
			sum[i] += 1- s/(double) K ;
		}
	}
	for (int i = 0; i< 16; i++)
		cout << sum[i] / Run << '\t';
	cout << '\n';
	
	int fTime = time(0);
	cout << "Time: " << fTime - startTiem<< '\n';
	
	r.RandomInit(seed);
	for (int i = 0; i< 16; i++)
		sum[i] = 0;
	
	
	for (int run=0; run<Run; run++) {
		//cout << "Run: " << run << "\n";
		LT_sim<Bit> sim;
		sim.init(K, K*1.2, 10, Degree_of_Edge, Omega, r.BRandom());
		
		
		
		for (int i = 0; i< 16; i++) {
			//			double t = Encoder(K, K*(1.05+0.01*i), Dsize);
			//			ErrorCount[i][(int)(t*16)]++;
			//			BER[i]+=t;
			
			
			
						sim.seqReceive(K*(1.05+0.01*i) -1);
						//sim.decode();
			
			
			
			//			for (int i=0; i< c.size(); i++) {
			//				if (r.Random() < 0.1) {
			//					c[i].setErased(1);
			//				}
			//			}
			
			//for (int j=K*(1.04+0.01*i); j<K*(1.05+0.01*i); j++) {
//				c[j].setErased(false);
//			}
//			Codeword<Bit> out = sim.decode(c);
//			
			
			sum[i] += sim.getNumOfErased()/(double)K;
			//cout << sim.run();
//			int s=0;
//			for (int j = 0; j<K; j++) {
//				if ( !out[j].isErased() && in[j] == out[j]) 
//				{
//					s++;
//				}
//				
//			}
//			sum[i] += 1- s/(double) K ;
		}
	}
	for (int i = 0; i< 16; i++)
		cout << sum[i] / Run << '\t';

	cout << "Time: " << time(0) - fTime << '\n';	
    return 0;
}
