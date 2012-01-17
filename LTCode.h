/*
 *  LTCode.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/2.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *	
 *	for a input codeword of length 3000000 would cost about 
 *	500MB of memory
 */

//#include "CodeSim.h"
#include "randomc.h"
#include "LT.h"
#include <limits>
//#include <hash_set.h>

#ifndef H_LTCode
#define H_LTCode

namespace CodeSim {
	
	/**
	 *	@brief Multiple-block LT coder
	 */
	template<class S>
	class LTCode : public CodingBlock<S,S> {
	public:
		LTCode(){}
		LTCode(int k, int max_n, int tag_size, int *tags, double * omega, int seed){
			this->k = k;
			this->max_n = max_n;
			this->tag_size = tag_size;
			this->tags = tags;
			this->omega = omega;
			this->seed = seed;
			reset();
		}
		
		void reset(){
			coders.clear();
		}
		Codeword<S> encode(Codeword<S>& a);
		Codeword<S> decode(Codeword<S>& a);
	private:
		
		int k, max_n, tag_size, *tags;
		double *omega;
		int seed;
		vector< LT_sim<S> > coders;
	};
	
	template<class S>
	Codeword<S> LTCode<S>::encode(Codeword<S>& a){
		Codeword<S> output;
		CRandomMersenne Rnd(seed);
		reset();
		for (int i=0; i<a.size(); i+=k) {
			coders.push_back(LT_sim<S>(k,max_n, tag_size, tags, omega, Rnd.BRandom()));
			Codeword<S> temp;
			if (i+k < a.size()) {
				temp.assign(a.begin()+i, a.begin()+i+k);
			}
			else{
				temp.assign(a.begin()+i, a.end());
				temp.insert(temp.end(), k-temp.size(), 0);
			}
			
			temp = coders.back().encode(temp);
			output.insert(output.end(), temp.begin(), temp.end());
		}
		
		stack<int> s = a.getMessageStack();
		s.push(a.size());
		s.push(seed);
		output.setMessageStack(s);
		return output;
	}
	
	template<class S>
	Codeword<S> LTCode<S>::decode(Codeword<S>& a){
		Codeword<S> output;
		//reset();
		stack<int> s = a.getMessageStack();
		if (s.top() != seed) {
			seed = s.top();
			reset();
		}
		CRandomMersenne Rnd(seed);
		for (int i=0; i<coders.size(); i++) {
			Rnd.BRandom();
		}
		for (int i=0; i<a.size(); i+=max_n) {
			
			Codeword<S> temp;
			if (i+max_n < a.size()) {
				temp.assign(a.begin()+i, a.begin()+i+max_n);
			}
			else{
				temp.assign(a.begin()+i, a.end());
				temp.insert(temp.end(), max_n-temp.size(), 0);
			}
			while(i/max_n > coders.size()) {
				coders.push_back(LT_sim<S>(k,max_n, tag_size, tags, omega, Rnd.BRandom()));
			}
			temp = coders[i/max_n].decode(temp);
			output.insert(output.end(), temp.begin(), temp.end());
		}
		
		s.pop();//pop seed
		output.trim(s.top());
		s.pop();
		output.setMessageStack(s);
		return output;
	}
	
}

#endif

