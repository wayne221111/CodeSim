/*
 *  LT.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <queue>
#include "randomc.h"
#include "CodeSim.h"
using namespace std;

#ifndef H_LT
#define H_LT

namespace CodeSim{

	/**
	 *	@brief Single block LT Coder
	 */
	template<class S>
	class LT_sim : public CodingBlock<S,S>{
	public:
		LT_sim();
		LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void init(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void seqReceive(int t);
		void receive(int t);
		
		bool isDecoded(int t);
		double failureRate();
		Codeword<S> encode(Codeword<S>& a);
		Codeword<S> decode(Codeword<S>& a);
		void reset();
		Codeword<S> getResult();
		void codeGen(int t);
		void printGraph(ostream &os);
		void printDecodingSequence(ostream &os);
		int getNumOfErased();
		
	private:
		
		void receive(int t, S s);
		void decode();
		int Num_of_Input, Num_of_Output, Num_of_Degree, Num_of_Decoding, 
				//ReceivedSize,
				generatedCode;
		int seed;
		CRandomMersenne Rnd;
		vector<long int>	d;//(Num_of_Output, 0);
		vector<long int>	output_symbol_degree;
		vector<long int>	decoding_degree_sequence;
		vector<vector<long int> >	edge;// = new vector<long int>[Num_of_Output];
		//vector<long int>	R_M;//(Num_of_Output, 0);
		//vector<long int>	erasure;//(Num_of_Output, 0);
		vector<char>	DE;//(Num_of_Input, 0);
		vector<int>			Degree_of_Edge;
		vector<double>		Omega;
		vector<int>			receivedMask;
		Codeword<S>			Result, Mid_Output;
		queue<long int>		ripple;
		vector< queue<long int> > input_edge; 
	};
	
	template<class S>
	LT_sim<S>::LT_sim():Rnd(0){}
	
	template<class S>
	LT_sim<S>::LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed):Rnd(seed){
		init(k, max_n, tag_size, tags, omega, seed);
	}
	template<class S>
	void LT_sim<S>::init(int k, int max_n, int tag_size, int *tags, double * omega, int seed){
		Num_of_Input = k;
		Num_of_Output= max_n;
		Num_of_Degree = tag_size;
		
		this->seed = seed;
		Degree_of_Edge.assign(tags, tags+tag_size);
		Omega.assign(omega, omega+tag_size);
		for(int i=1; i<Num_of_Degree; i++)  
			Omega[i] = Omega[i-1]+Omega[i];
		//Normalize
		for(int i=1; i<Num_of_Degree; i++)  
			Omega[i] /= Omega[Num_of_Degree-1];
		
		reset();

	}
	
	template<class S>
	bool LT_sim<S>::isDecoded(int t){
		return DE[t];
	}
	
	template<class S>
	void LT_sim<S>::codeGen(int t){
		if (generatedCode < t+1){
			for(int i=generatedCode; i<=t; i++){
				// decide degree
				double rando = Rnd.Random();
				int s, flag;
				
				if(rando<Omega[0])
				{
					d[i] = Degree_of_Edge[0];
				}      
				
				else
				{
					for(s=0;s<Num_of_Degree-1;s++)
					{
						if(rando>=Omega[s] && rando<Omega[s+1])
						{
							d[i] = Degree_of_Edge[s+1];
							break;
						}
					}
				}
				output_symbol_degree[i] = d[i];
				if (d[i] == 0) {
					cerr << "Error: node of degree 0.\n";
				}
				// decide connection
				edge[i].assign( d[i], 0 );
				for(int j=0;j<d[i];j++)
				{
					flag=1;
					while(flag==1)
					{
						flag=0;
						//	edge[i][j]=rand()%Num_of_Input;
						edge[i][j]=Rnd.IRandomX(0,Num_of_Input-1);
						for(int lllll=0;lllll<j;lllll++)
						{
							if(edge[i][j]==edge[i][lllll])
							{
								flag=1;
							}
						}
					}
				}
				
				
			}
			
			generatedCode = t+1;
		}
	}
	
	template<class S>
	inline void LT_sim<S>::receive(int t){
		receive(t, 0);
	}
	
	template<class S>
	void LT_sim<S>::seqReceive(int t){
		for(int i=0; i <= t; i++)
		{
			receive(i);
		}
	}
	
	template<class S>
	inline void LT_sim<S>::receive(int t, S s){
		if (Num_of_Decoding >= Num_of_Input) {
			return;
		}
		codeGen(t);
		if(receivedMask[t]==0)
		{
			receivedMask[t]=1;
			for (int i=0; i<d[t]; i++) {
				if (DE[ edge[t][i] ] == 1) {
					s = s + Result[ edge[t][i] ];
					d[t]--;
					edge[t][i] = edge[t][d[t]];
					i--;
				}
				else {
					input_edge[ edge[t][i] ].push(t);
				}

			}
			
			if (d[t] == 0) {
				return;
			}
			if (d[t] == 1) {
				decoding_degree_sequence.push_back(output_symbol_degree[t]);
				ripple.push(edge[t][0]);
				DE[ edge[t][0] ] = 1;
				Result[ edge[t][0] ] = s;
				Num_of_Decoding++;
				return;
			}
			
			
			//R_M[ReceivedSize]=t;
			//ReceivedSize++;
			
			Mid_Output[t] = s;
		}
		
	}
	
	template<class S>
	void LT_sim<S>::decode(){
		//int flag=1;
		while( Num_of_Decoding < Num_of_Input && !ripple.empty())
		{
			long int r = ripple.front();
			
			while (!input_edge[r].empty()) {
				long int t = input_edge[r].front();
				for (int i=0; i<d[t]; i++) {
					if(edge[t][i] == r){
						Mid_Output[t] = Mid_Output[t] + Result[r];
						d[t]--;
						edge[t][i] = edge[t][d[t]];
						break;
					}
				}
				if (d[t] == 1 && DE[ edge[t][0] ] == 0) {
					decoding_degree_sequence.push_back(output_symbol_degree[t]);
					DE[ edge[t][0] ]=1;
					ripple.push(edge[t][0]);
					Result[ edge[t][0] ] = Mid_Output[t];
					Num_of_Decoding++;
				}
				input_edge[r].pop();
			}
			
			
			ripple.pop();
//			flag=0;
//			for(int i=0;i<ReceivedSize;i++)
//			{
//				if(d[R_M[i]]==1 && DE[ edge[ R_M[i] ][0] ]==0)
//				{
//					DE[ edge[ R_M[i] ][0] ]=1;
//					Result[ edge[R_M[i]][0] ] = Mid_Output[R_M[i]];
//					Result[ edge[ R_M[i] ][0] ].setErased(false);
//					Num_of_Decoding++;
//					ReceivedSize--;
//					R_M[i]=R_M[ReceivedSize];
//					i--;
//					flag=1;
//				}
//			}
//			
//			for(int i=0;i<ReceivedSize;i++)
//			{
//				for(int j=0;j<d[R_M[i]];j++)
//				{
//					if(DE[ edge[ R_M[i] ][j] ]==1)
//					{
//						Mid_Output[R_M[i]] = Mid_Output[R_M[i]] + Result[edge[R_M[i]][j]];
//						d[R_M[i]]--;
//						edge[ R_M[i] ][j]=edge[ R_M[i] ][ d[R_M[i]] ];
//						j--;
//					}
//				}
//			}
		} 
	}
	
	template<class S>
	Codeword<S> LT_sim<S>::encode(Codeword<S>& a){
		reset();
		codeGen(Num_of_Output-1);
		stack<int> s = a.getMessageStack();
		Codeword<S> output;
		output.assign(Num_of_Output,0);
		
		for (int i=0; i< Num_of_Output; i++) {
			for (int j=0; j<edge[i].size(); j++) {
				output[i] = output[i] + a[edge[i][j]];
			}
		}
		
		Mid_Output = output;
		s.push(seed);
		output.setMessageStack(s);
		
		return output;
	}
	
	template<class S>
	Codeword<S> LT_sim<S>::decode(Codeword<S>& a){
		stack<int> s = a.getMessageStack();
		if (s.size() != 0 && s.top() != seed) {
			seed = s.top();
			reset();
			codeGen(Num_of_Output-1);
			Mid_Output = a;
		}
		
//		if (ReceivedSize == 0) {
//			Mid_Output = a;
//		}
		
		for (int i =0; i< a.size() && i< Num_of_Output; i++) {
			if (!a[i].isErased() ) {
				receive(i, a[i]);
			}
		}
		
		
		decode();
		if (s.size()>0) {
			s.pop();
			Result.setMessageStack(s);
		}
		return Result;
	}
	
	template<class S>
	void LT_sim<S>::reset(){
		//	Num_of_Input = k;
		//	Num_of_Output= max_n;
		//	Num_of_Degree = tag_size;
		
		d.assign(Num_of_Output,0);
		output_symbol_degree.assign(Num_of_Output,0);
		decoding_degree_sequence.clear();
		edge.assign(Num_of_Output,vector<long int>());
		//R_M.assign(Num_of_Output,0);
		//erasure.assign(max_n,0);
		
		DE.assign(Num_of_Input,0);
		Rnd.RandomInit(seed);
		
		
		Num_of_Decoding = 0;
		//ReceivedSize = 0;
		generatedCode = 0;
		//codeGen(Num_of_Output);
		
		receivedMask.assign(Num_of_Output, 0);
		
		Result.assign(Num_of_Input, -1);
		Mid_Output.assign(Num_of_Output, 0);
		
		ripple = queue<long int>();
		input_edge.assign(Num_of_Input, queue<long int>());
	}
	
	template<class S>
	double LT_sim<S>::failureRate(){
		decode();
		return 1-(Num_of_Decoding/(double)Num_of_Input);//failure rate
		
	}
	template<class S>
	Codeword<S> LT_sim<S>::getResult(){
		decode();
		return Result;
	}
	
	template<class S>
	void LT_sim<S>::printGraph(ostream &os){
		os << "*Nodes\nid*int	label*string	symbol_type*string\n";
		for (int i=0; i<Num_of_Input; i++) {
			os << i+1 << "\t\"source_" << i+1 << "\"\t\"s\""  << '\n';
		}
		
		for (int i=0; i<generatedCode; i++) {
			os << Num_of_Input+i+1 << "\t\"encoding_" << i+1 << "\"\t\"e\"" << '\n';
		}
		
		os << "*UndirectedEdges\nsource*int	target*int\n";
		
		for (int i=0; i<generatedCode; i++) {
			for (int j=0; j<d[i]; j++) {
				os << Num_of_Input+i+1 << '\t' << edge[i][j] +1 << '\n';
			}
		}
	}
	
	template<class S>
	void LT_sim<S>::printDecodingSequence(ostream &os){
		for (int i=0; i<decoding_degree_sequence.size(); i++) {
			os << decoding_degree_sequence[i] << '\t';
		}
		os << endl;
	}
	
	template<class S>
	int LT_sim<S>::getNumOfErased() {
		decode();
		return Num_of_Input - Num_of_Decoding;
	}
}

#endif
