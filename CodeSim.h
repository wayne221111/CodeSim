/*
 *  CodeSim.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/25.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include<vector>
#include<string>
#include<stack>
#include<fstream>
#include<iostream>
#include<cstdlib>
using namespace std;

#ifndef H_CodeSim_BASIC
#define H_CodeSim_BASIC

#define PARALLEL_THREADS 10

namespace CodeSim {
	
	template<class T, int W>
	class Symbol{
	public:
		Symbol(){
			erased = false;
		}
		
		
		bool isErased(){
			return erased;
		}
		
		void setErased(bool t){
			erased = t;
		}
		
		T getValue(){
			return value;
		}
		
		bool operator==(Symbol t){
			return this->value == t.value;
		}
		
	protected:
		T value:W;
	private:
		bool erased:1;
	};
	
	
	
	class Bit : public Symbol<bool,1>{
	public:
		Bit();
		
		
		Bit(int t);
		

		Bit operator+(Bit t);		
		//Bit operator*(Bit t);		
		string toString();
	};
	
	class Byte : public Symbol<unsigned char,8>{
	public:
		static const unsigned int mask = (1<<8)-1;
		Byte();
		
		
		Byte(int t);
		
		
		Byte operator+(Byte t);		
		//Byte operator*(Byte t);		
		string toString();
	};
	
	template<class S>
	class Codeword : public vector<S> {
	public:
		Codeword(){};
		Codeword(int size, S s){
			assign(size,s);
		}
		stack<int>&  getMessageStack(){
			return messageStack;
		}
		
		void setMessageStack(stack<int> s){
			messageStack = s;
		}
		
		void trim(int m){
			if(m < this->size())
			{
				this->erase(this->begin()+m, this->end());
			}
		}
	private:
		stack<int> messageStack;
	};
	
	template<class S1, class S2>
	class CodingBlock {
	public:
		
//		virtual Codeword<S2> forward(Codeword<S1> c) = 0;
//		virtual Codeword<S1> backward(Codeword<S2> c) = 0;
		
	private:
	};
	
	class BitToByteCoverter : public CodingBlock<Bit, Byte> {
	public:
		 static Codeword<Byte> convert(Codeword<Bit> &a){
			Codeword<Byte> output;
			output.reserve(a.size()/8+1);
			stack<int> s = a.getMessageStack();
			s.push(a.size());
			
			int i;
			for (i=0; i<a.size(); i+=8) {
				int v=0;
				bool e=false;
				
				for (int j=0; j<8; j++) {
					v <<= 1;
					if(i+j >= a.size())
						continue;
					
					if(a[i+j].getValue())
						v ^= 1;
					if(a[i+j].isErased())
						e=true;
				}
				
				if(e){
					output.push_back(-1);
				}
				else {
					output.push_back(v);
				}

			}
			
			output.setMessageStack(s);
			return output;
		}
		
		static Codeword<Bit> revert(Codeword<Byte> &a){
			Codeword<Bit> output;
			output.reserve(a.size()*8);
			stack<int> s = a.getMessageStack();
			
			for (int i=0; i<a.size(); i++) {
				if (a[i].isErased()) { // Byte is erased
					// 8 erased Bits
					for (int i=0; i<8; i++) {
						output.push_back(-1); 
					}
				}
				else {
					int array[8];
					unsigned int tmp = a[i].getValue();
					for (int i=0; i<8; i++) {
						if (tmp & 1){
							array[7-i] = 1;
						}
						else {
							array[7-i] = 0;
						}

						tmp >>= 1;
					}
					
					output.insert(output.end(), array, array+8);
				}

			}
			
			output.trim(s.top());
			s.pop();
			output.setMessageStack(s);
			
			return output;
		}
	};
	
	template<class T>
	class Permutator : public CodingBlock<T,T> {
	public:
		
		// default construct, blockSize = 1 (no permutation).
		Permutator(){
			blockSize = 1;
			permutationTable.assign(1,0);
			depermutationTable.assign(1,0);
		}
		
		/* 
			read permutation table from file "filename"
			inverseOrder = true would inverse the table (permutate() become depermutate(), and vice versa)
		*/
		void init(const string filename, bool inverseOrder){
			ifstream in(filename.c_str());
			if(in.fail())
			{
				cerr << "Permutator initial error: \""<< filename << "\" can't be opened." ;
				exit(-1);
				return;
			}
			
			
			depermutationTable.clear();
			
			
			permutationTable.clear();
			unsigned int t;
			while (in >> t) {
				permutationTable.push_back(t);
			}
			
			blockSize = permutationTable.size();
			// validate
			vector<bool> checked(blockSize, 0);
			for (int i = 0; i< permutationTable.size(); i++) {
				t = permutationTable[i];
				if( t < blockSize) {
					checked[t] = 1;
				}
				else {
					cerr << "Permutator initial error: \""<< filename << "\" is not a valid interleaver: index out of range."<<endl ;
					exit(-1);
					return;
				}

			}
			
			for (int i = 0; i< checked.size(); i++) {
				if(checked[i] == 0)
				{
					cerr << "Permutator initial error: \""<< filename << "\" is not a valid interleaver: duplicated indexes."<<endl ;
					exit(-1);
					return;
				}
			}
			
			depermutationTable.assign(blockSize, 0);
			
			for (int i = 0; i< permutationTable.size(); i++) {
				t = permutationTable[i];
				depermutationTable[t] = i;
			}
			
			
			if(inverseOrder){
				permutationTable.swap(depermutationTable);
			}
			
			
		}
		
		/** 
		 * @brief read permutation table from file "filename"
		 * @param
		 */
		
		Permutator(const string filename){
			init(filename, false);
		}
		
		Permutator(const string filename, bool inverseOrder){
			init(filename, inverseOrder);
		}
		
		Codeword<T> permutate(Codeword<T> c){
			Codeword<T> output;
			for (int i = 0; i< c.size(); i+=blockSize) {
				Codeword<T> temp;
				temp.assign(blockSize, 0);
				for (int j=0; j< blockSize && i+j < c.size(); j++) {
					temp[permutationTable[j]] = c[i+j];
				}
				
				output.insert(output.end(), temp.begin(), temp.end());
			}
			stack<int> s = c.getMessageStack();
			s.push(c.size());
			output.setMessageStack(s);
			return output;
		}
		
		Codeword<T> depermutate(Codeword<T> c){
			Codeword<T> output;
			for (int i = 0; i< c.size(); i+=blockSize) {
				Codeword<T> temp;
				temp.assign(blockSize, 0);
				for (int j=0; j< blockSize && i+j < c.size(); j++) {
					temp[depermutationTable[j]] = c[i+j];
				}
				
				output.insert(output.end(), temp.begin(), temp.end());
			}
			stack<int> s = c.getMessageStack();
			if (s.size()>0) {
				output.trim(s.top());
				s.pop();
				output.setMessageStack(s);
			}
			return output;
		}
		
		Codeword<T> forward(Codeword<T> c){
			return permutate(c);
		}
		
		Codeword<T> backward(Codeword<T> c){
			return depermutate(c);
		}
		
		
		void print(){
			for(int i = 0; i< blockSize; i++){
				cout << depermutationTable[i] << endl;
			}
		}
	private:
		int blockSize;
		vector<unsigned int> permutationTable, depermutationTable;
	};
	
	
	
}
#endif

