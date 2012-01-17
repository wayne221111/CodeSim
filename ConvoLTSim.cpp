/*
 *  ConvoLTSim.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/5/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "ConvoCode.h"
#include "LTCode.h"
#include "randomc.h"
#include <ctime>
#include <cmath>
#include <map>
//#define L 100000
#define MAX_BIT 1000000000ULL
#define BASE 1.08
#define STEPS 1
#define Delta 0.01
#include <omp.h>

using namespace CodeSim;
using namespace std;

int main(int argn, char **args){
	
	if (argn < 4) {
		cerr << "Usage: ConvoLTSim.out convo_code_file LT_file interleaver_file" << endl;
		exit(-1);
	}
	
	int start_time = time(0);
//	string puncher_table = "1";
//	
//	
//	if (argn > 3) {
//		puncher_table = args[3];
//	}
	
	CRandomMersenne r(time(0));
	ConvoCode cc( args[1] );
	int Layer = cc.getK();
	unsigned long L = 80000*6/7/Layer;
	unsigned long Run = MAX_BIT / L;
	ifstream ifsLT(args[2]);
	if(ifsLT.fail()){
		cerr << "Error: file \""+string(args[2])+"\" does not exist." << endl;
		exit(-1);
	}
	int K, MaxN, Dsize, *Tags;
	double *Distribution;
	
	ifsLT >> K >> Dsize;
	
	Tags = new int[Dsize];
	Distribution = new double[Dsize];
	
	for (int i=0; i<Dsize; i++) {
		ifsLT >> Tags[i];
	}
	for (int i=0; i<Dsize; i++) {
		ifsLT >> Distribution[i];
	}
	
	
	Permutator<Bit> inter(args[3], true);
	
	
	cout << "ConvoCode file: \"" << args[1] <<"\" LT file: " << args[2] << "\" interleaver file: " << args[3] << '"' << endl;
	cout << "Epsilon";
	for (int i=0; i<Layer; i++) {
		cout << "\tLayer" << i+1;
	}
	cout << "\tTotalBits\tLT\tLT_total\tRun\n";
	
	
	for (int s=0; s<STEPS; s++) {
		vector<unsigned long> total_err(Layer, 0);
		unsigned long total=0, LT_total = 0, LT_total_err = 0;
		double Epsilon = (BASE+s*Delta) -1;
		vector<vector<bool> > histo_mask(Layer, vector<bool>(L+1,0)), LT_histo_mask(1, vector<bool>(L+1,0));
		vector<map<unsigned long, unsigned long> > histo(Layer, map<unsigned long, unsigned long>()), LT_histo(1, map<unsigned long, unsigned long>());
		//		for (int i=0; i<Layer; i++) {
		//			err[i]=0;
		//		}
		
		//while (1) 
		{
			#pragma omp parallel for schedule(dynamic) num_threads(6)
			for (int i=0; i<Run; i++) {
				vector<unsigned long> err(Layer, 0);
				unsigned long LT_err = 0;
				Codeword<Bit> a;
				a.reserve(Layer*L);
				for (int t=0; t<Layer*L; t++) {
					a.push_back(r.IRandomX(0, 1));
				}
				Codeword<Bit> b = cc.encode(a);
				b = inter.permutate(b);
				Codeword<Byte> b1 = BitToByteCoverter::convert(b);
				
				if(s == 0 && i==0)
					cout << "LT block:" << b1.size() << endl;
				LT_sim<Byte> lt(b1.size(), b1.size()*(1+Epsilon), Dsize, Tags, Distribution, r.BRandom());
				Codeword<Byte> c1 = lt.encode(b1);
				c1 = lt.decode(c1);
				for (int i=0; i<c1.size(); i++) {
					if(c1[i].isErased()) {
						LT_err++;
					}
				}
				
				#pragma omp atomic
				LT_total_err+= LT_err;
				#pragma omp atomic
				LT_total += c1.size();
				#pragma omp critical
				{
					if(LT_histo_mask[0][LT_err]){
						LT_histo[0][LT_err] ++;
					}
					else {
						LT_histo[0][LT_err] = 1;
						LT_histo_mask[0][LT_err] =1;
					}
				}
				
				if (LT_err > 0)
				{
					Codeword<Bit> c2 = BitToByteCoverter::revert(c1);
					c2 = inter.depermutate(c2);

					Codeword<Bit> c = cc.decode(c2);

					for (int i=0; i<c.size(); i++) {
						if (!(a[i] == c[i])) {
							//#pragma omp atomic
							err[i%Layer]++;
						}
					}
				}
				
				#pragma omp atomic
				total += L;
				
				for (int i=0; i<Layer; i++) {
					#pragma omp atomic
					total_err[i]+=err[i];
				}
				
				#pragma omp critical
				for (int i=0; i<Layer; i++) {
					if(histo_mask[i][err[i]]){
						histo[i][err[i]] ++;
					}
					else {
						histo[i][err[i]] = 1;
						histo_mask[i][err[i]] =1;
					}
					
				}
			}
			
			
		}
		
		// enough # of errors
		
		cout << Epsilon;
		
		for (int i=0; i<Layer; i++) {
			cout << '\t' << total_err[i];
		}
		
		cout << '\t' << total << '\t' << LT_total_err << '\t' << LT_total << '\t' << Run << endl;
		
		for (int i=0; i<Layer; i++) {
			cout << "Layer " << i +1 << " histogram:\n";
			for (map<unsigned long, unsigned long>::iterator it = histo[i].begin(); it!=histo[i].end(); it++) {
				cout << it->first << '\t' << it->second  << '\n';
			}
		}
		
		cout << endl;
		
		cout << "Lt histogram:\n";
		for (map<unsigned long, unsigned long>::iterator it = LT_histo[0].begin(); it!=LT_histo[0].end(); it++) {
			cout << it->first << '\t' << it->second << '\n';
		}
	}
	
	cout << "Time: " << time(0)-start_time << endl;
	
}