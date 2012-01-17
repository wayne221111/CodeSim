/*
 *  ConvoCodeSim_byte.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/6/1.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "ConvoCode.h"
#include "LTCode.h"
#include "randomc.h"
#include <ctime>
#include <cmath>
#define L 100000
#define MAX_BIT 1000000000L
#define MIN_BIT 10000000L
#define START_ERROR_RATE 0.1
#define STEPS_PER_ORDER 1
#define MIN_ERRORS 100
#include <omp.h>

using namespace CodeSim;
using namespace std;

int main(int argn, char **args){
	
	if (argn <= 2) {
		cerr << "Usage: ConvoCodeSim_byte.out code_file interleaver_file" << endl;
		exit(-1);
	}
	
	
	//string puncher_table = "1";
	
	
//	if (argn > 2) {
//		puncher_table = args[2];
//	}
	
	CRandomMersenne r(time(0));
	ConvoCode cc( args[1] );
	Permutator<Bit> inter(args[2], true);
	int Layer = cc.getK();
	cout << "Code file: \"" << args[1] << "\" Interleaver file: \""<< args[2] << '"' << endl;
	cout << "ChannelErasureRate";
	for (int i=0; i<Layer; i++) {
		cout << "\tLayer" << i+1;
	}
	cout << "\tTotalBits\n";
	
	double errRate = START_ERROR_RATE;
	while (1) {
		vector<unsigned long> err(Layer, 0);
		unsigned long total=0;
		
		//		for (int i=0; i<Layer; i++) {
		//			err[i]=0;
		//		}
		
		while (1) {
			#pragma omp parallel for num_threads(6)
			for (int i=0; i<6; i++) {
				Codeword<Bit> a;
				a.reserve(Layer*L);
				for (int i=0; i<Layer*L; i++) {
					a.push_back(r.IRandomX(0, 1));
				}
				Codeword<Bit> a2 = cc.encode(a);
				a2 = inter.permutate(a2);
				Codeword<Byte> b = BitToByteCoverter::convert(a2);
				
				for (int i=0; i< b.size(); i++) { 
					if (r.Random() < errRate) {
						b[i].setErased(true);
					}
				}
				
				Codeword<Bit> b2 = BitToByteCoverter::revert(b);
				b2 = inter.depermutate(b2);
				Codeword<Bit> c = cc.decode(b2);
				
				for (int i=0; i<c.size(); i++) {
					if (!(a[i] == c[i])) {
						#pragma omp atomic
						err[i%Layer]++;
					}
				}
				#pragma omp atomic
				total += L;
			}
			// check stop
			
			if (total < MIN_BIT) {
				continue;
			}
			if (total >= MAX_BIT) {
				break;
			}
			
			bool enough = true;
			for (int i=0; i<Layer; i++) {
				if (err[i] < MIN_ERRORS) {
					enough = false;
				}
			}
			
			if (enough) {
				break;
			}
		}
		
		// enough # of errors
		
		cout << errRate;
		
		for (int i=0; i<Layer; i++) {
			cout << '\t' << err[i]/(double)total;
		}
		
		cout << '\t' << total << endl;
		
		errRate *= pow(10, -1.0/STEPS_PER_ORDER);
		
	}
	
}
