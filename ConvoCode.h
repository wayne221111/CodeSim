/*
 *  ConvoCode.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "CodeSim.h"
using namespace std;

#ifndef H_ConvoCode
#define H_Convocode


namespace CodeSim {
	int octToDec(int t);
	int setBit(int i, int position, bool value);
	unsigned int countBit(unsigned int);
	
	
	template<class T>
	void intToArray(T *array, unsigned int value, int width){
		for(int j=0; j< width; j++){
			array[width-1-j] = (value&1);
			value >>= 1;
		}
	}
	
	
	struct state {
		int out;
		int next;
		state(){
			out = 0;
			next = 0;
		}
		state(int output, int nextState){
			out = output;
			next = nextState;
		}
	};
	
	
	class ConvoCode : public CodingBlock<Bit,Bit> {
	public:
		ConvoCode();
		ConvoCode(string filename);
		void showInfo();
		Codeword<Bit> encode(Codeword<Bit> a) const;
		Codeword<Bit> decode(Codeword<Bit> & a) const;
		int getK() const;
		int getN() const;
	protected:
		void generateTrellis();
	private:
		int k, n, m, max_forney;
		vector<int> forneyIndices;
		vector< vector<int> > G;
		vector< vector<state> > trellis;
		string puntureTable;
	};

}
#endif
