/*
 *  ConvoCode.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/7.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "ConvoCode.h"
#include <list>

namespace CodeSim {

	int octToDec(int t){
		if(t <= 0)
			return 0;
		else {
			return ((t%10)&7) ^ (octToDec(t/10) << 3);
		}
		
	}
	
	int setBit(int i, int position, bool value){
		i = i & (-1 - (1 << position) );
		int t = value;
		t <<= position;
		return i+t;
	}
	
	unsigned int countBit(unsigned int i){
		if(i == 0)
			return 0;
		if(i&1)
			return 1 + countBit(i>>1);
		else {
			return countBit(i>>1);
		}
		
	}
	
	
	
	/**
	 *	@brief Construct ConvoCode from file.
	 *	@param filename The file to be read.
	 *
	 *
	 *	The file should consist 2+k lines.
	 *	1st line: N and K
	 *	2nd line: Forney Indices(K numbers)
	 *	3rd to 2+K line: generator matrix (K lines of N numbers, int octal)
	 */
	
	ConvoCode::ConvoCode(string filename){
		ifstream fin(filename.c_str());
		if(fin.fail()){
			cerr << "File opening Error: " << filename << '\n';
		}
		
		fin>>n>>k;
		forneyIndices.assign(k, 0 );
		m=0;
		max_forney = 0;
		for (int i=0; i<k; i++) {
			fin >> forneyIndices[i];
			m += forneyIndices[i];
			if(forneyIndices[i] > max_forney)
				max_forney = forneyIndices[i];
		}
		
		G.assign(k,vector<int>() );
		for (int i=0; i<k; i++) {
			G[i].assign(n,0);
			for(int j=0; j<n; j++){
				int t;
				fin >> t;
				G[i][j] = octToDec(t);
			}
		}
		
		if(!(fin >> puntureTable))
			puntureTable = "1";
		fin.close();
		generateTrellis();
		
	}
	
	void ConvoCode::showInfo(){
		cout << "N: " << n << "\tK: " << k 
		<< "\nForney Indices: ";
		for (int i = 0; i< forneyIndices.size(); i++)
			cout << '\t' <<forneyIndices[i];
		cout << "\nG:\n";
		for(int i = 0; i< k; i++){
			for (int j=0; j< n; j++){
				cout << G[i][j] << ' ';
			}
			cout << '\n';
		}
		cout << "Trellis:\n";
		for(int i=0; i< trellis.size(); i++)
		{
			for(int j=0;j<trellis[i].size(); j++){
				cout << trellis[i][j].out << ' ';
			}
			cout << '\n';
		}
		
	}
	
	void ConvoCode::generateTrellis(){
		trellis.assign(1<<m, vector<state>());
		for(int i=0; i< (1<<m); i++){
			trellis[i].assign(1<<k, state());
			for(int j=0; j< (1<<k); j++){
				// determine next state
				int t = i << 1, p=0, input = j; ;
				for(int d=0; d<k; d++){
					if(forneyIndices[k-1-d] > 0){ 
						t = setBit(t, p, input & 1);
						p += forneyIndices[k-1-d];
					}
					input >>= 1;
				}
				trellis[i][j].next = t & ((1<<m) -1);
				
				// determine output
				t = 0, p=0, input = j;
				for(int d=0; d<k; d++){
					if(forneyIndices[k-1-d] > 0){ 
						t += (( i >> (p-d) ) & ( (1<<forneyIndices[k-1-d]) -1)) << (p+1); // take out the memory of (k-1-d)th input
					}
					t = setBit(t, p, input & 1);
					p += forneyIndices[k-1-d]+1;
					input >>= 1;
				}
				
				
				for(int column = 0; column < n; column++){
					// get column of G
					int g = 0;
					trellis[i][j].out <<= 1;
					for(int row=0; row < k; row++){
						g <<= forneyIndices[row]+1;
						g += G[row][column];
					}
					trellis[i][j].out += countBit(g & t) & 1;
				}
			}
			
		}
	}
	
	
	
	Codeword<Bit> ConvoCode::encode(Codeword<Bit> a) const{
		Codeword<Bit> output;
		output.reserve( (a.size()/k+1) *n);
		stack<int> s = a.getMessageStack();
		s.push(a.size());
		
		int i, *o = new int[n], m = 0;
		
		for(i=0; i+k<a.size(); i+=k){
			int t=0;
			for(int j=0; j<k; j++){
				t <<= 1;
				t += a[i+j].getValue();
			}
			
			intToArray(o, trellis[m][t].out, n);
			
			output.insert(output.end(), o, o+n);
			
			m = trellis[m][t].next;
			
		}
		// last bits
		int t=0;
		for(int j=0; j< k; j++){
			t <<= 1;
			if(i+j < a.size())
			{
				t += a[i+j].getValue();
			}
			
		}
		
		intToArray(o, trellis[m][t].out, n);
		
		output.insert(output.end(), o, o+n);
		
		m = trellis[m][t].next;
		// end last bits
		
		// todo closing
		for(i = 0; i< max_forney; i++){
			int t=0;
						
			intToArray(o, trellis[m][t].out, n);
			
			output.insert(output.end(), o, o+n);
			
			m = trellis[m][t].next;
		}
		
		delete [] o;
		output.setMessageStack(s);
		
		double prate=0;
		for(int i=0; i<puntureTable.size(); i++){
			if (puntureTable[i] == '1') {
				prate++;
			}
		}
		prate /= puntureTable.size();
		Codeword<Bit> tmp;
		tmp.setMessageStack(s);
		tmp.reserve(output.size()*prate + n);
		//punture
		int j=0;
		Codeword<Bit>::iterator it=output.begin();
		while( it!=output.end() ) {
			if (puntureTable[j%puntureTable.size()] == '1') {
				tmp.push_back(*it);
			}
			
			it++;
			j++;
		}
		return tmp;
	}
	
	Codeword<Bit> ConvoCode::decode(Codeword<Bit> & pa) const{
		double prate=0;
		for(int i=0; i<puntureTable.size(); i++){
			if (puntureTable[i] == '1') {
				prate++;
			}
		}
		prate /= puntureTable.size();
		
		Codeword<Bit> output, a;
		output.reserve( (a.size()/n+1) *k);
		a.reserve(pa.size()/prate + n);
		stack<int> s = pa.getMessageStack();
		
		//depunture
		int j=0, it=0;
		while( it< pa.size() ) {
			while (puntureTable[j%puntureTable.size()] == '0') {
				a.push_back(-1);
				j++;
			}
			a.push_back(pa[it]);
			it++;
			j++;
		}
		
		while (a.size() % n != 0) {
			a.push_back(-1);
		}
		list< vector<unsigned int> > preState, cost, preIn;
		
		preState.push_back( vector<unsigned int>(1<<m, 0) );
		
		cost.push_back( vector<unsigned int>(1<<m, 999999) );
		preIn.push_back( vector<unsigned int>(1<<m, 0) );
		
		cost.front()[0] = 0;
		
		int *o = new int[n], trellisLength = 0;
		for(int i=0; i< a.size()-max_forney*n ; i+=n){
			preState.push_back( vector<unsigned int>(1<<m, 0) );
			cost.push_back( vector<unsigned int>(1<<m, 99999999) );
			preIn.push_back( vector<unsigned int>(1<<m, 0) );
			
			//int p = i/n;
			for(int fromState=0; fromState< (1 << m); fromState++){
				for(int in=0; in < (1 << k); in++){
					intToArray(o, trellis[fromState][in].out, n);
					int tempCost = 0;
					for(int s=0; s<n; s++){
						if (a[i+s].isErased()) {
							tempCost += 1;
						}
						else if(a[i+s].getValue() != o[s]){
							tempCost += 9999;
						}
					}
					list< vector<unsigned int> >::reverse_iterator iCost = cost.rbegin();
					iCost++;
					if ( ( *cost.rbegin() )[ trellis[fromState][in].next ] > ( *iCost )[fromState] + tempCost ) {
						( *cost.rbegin() )[ trellis[fromState][in].next ] = ( *iCost )[fromState] + tempCost;
						preState.back()[ trellis[fromState][in].next ] = fromState;
						preIn.back()[ trellis[fromState][in].next ] = in;
					}
				}
			}
			
			unsigned int t_min = -1;
			for (int j=0; j< (1<<m); j++) {
				if(( *cost.rbegin() )[j] < t_min){
					t_min = ( *cost.rbegin() )[j];
				}
			}
			for (int j=0; j< (1<<m); j++) {
				( *cost.rbegin() )[j] -= t_min;
			}
			
			trellisLength++;
			if (trellisLength > 80) {
				//sliding window decoding
				int state = 0, windowSize = 0;
				list<unsigned int> out_temp;
				for (list< vector<unsigned int> >::reverse_iterator 
					 iState = preState.rbegin(), 
					 iStatePre = preState.rbegin(), 
					 iIn = preIn.rbegin()
					 ; ++iStatePre != preState.rend(); iState++, iIn++) {
					
					windowSize++;
					if (windowSize >= 10*m) {
						out_temp.push_front((*iIn)[state]);
					}
					state = (*iState)[state];
					
					
				}
				
				
				for (; !out_temp.empty(); out_temp.pop_front()) {
					intToArray(o, out_temp.front(), k);
					output.insert(output.end(), o, o+k);
					
					preState.pop_front();
					cost.pop_front();
					preIn.pop_front();
					trellisLength--;
				}
			}
			
		}
		
		// ending stage
		
		for(int i=a.size()-max_forney*n; i< a.size(); i+=n){
			preState.push_back( vector<unsigned int>(1<<m, 0) );
			cost.push_back( vector<unsigned int>(1<<m, 99999999) );
			preIn.push_back( vector<unsigned int>(1<<m, 0) );
			
			
			for(int fromState=0; fromState< (1 << m); fromState++){
				int in=0; 
				{
					intToArray(o, trellis[fromState][in].out, n);
					int tempCost = 0;
					for(int s=0; s<n; s++){
						if (a[i+s].isErased()) {
							tempCost += 1;
						}
						else if(a[i+s].getValue() != o[s]){
							tempCost += 9999;
						}
					}
					list< vector<unsigned int> >::reverse_iterator iCost = cost.rbegin();
					iCost++;
					if ( ( *cost.rbegin() )[ trellis[fromState][in].next ] > ( *iCost )[fromState] + tempCost ) {
						( *cost.rbegin() )[ trellis[fromState][in].next ] = ( *iCost )[fromState] + tempCost;
						preState.back()[ trellis[fromState][in].next ] = fromState;
						preIn.back()[ trellis[fromState][in].next ] = in;
					}
				}
			}
			
			unsigned int t_min = -1;
			for (int j=0; j< (1<<m); j++) {
				if(( *cost.rbegin() )[j] < t_min){
					t_min = ( *cost.rbegin() )[j];
				}
			}
			for (int j=0; j< (1<<m); j++) {
				( *cost.rbegin() )[j] -= t_min;
			}
			
		}
		
		
		// go through trellis(last time)
		int state = 0;
		list<unsigned int> out_temp;
		for (list< vector<unsigned int> >::reverse_iterator 
				iState = preState.rbegin(), 
				iStatePre = preState.rbegin(), 
				iIn = preIn.rbegin()
				; ++iStatePre != preState.rend(); iState++, iIn++) {
			
			out_temp.push_front((*iIn)[state]);
			state = (*iState)[state];
		}
		
		
		for (; !out_temp.empty(); out_temp.pop_front()) {
			intToArray(o, out_temp.front(), k);
			output.insert(output.end(), o, o+k);
		}
		
		
		delete [] o;
		output.trim(s.top());
		s.pop();
		output.setMessageStack(s);
		return output;
	}
	
	int ConvoCode::getK() const{
		return k;
	}
	int ConvoCode::getN() const{
		return n;
	}
}
