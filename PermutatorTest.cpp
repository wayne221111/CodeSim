/*
 *  PermutatorTest.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/29.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include<iostream>
#include<string>
#include "CodeSim.h"
using namespace CodeSim;
using namespace std;
int main(){
	string file = "Interleaver2.txt";
	Permutator<int> p(file, false);
	p.print();
	
	Codeword<int> c;
	for (int i = 1; i<=10; i++) {
		c.push_back(i);
	}
	
	for (int i=0; i<c.size(); i++) {
		cout << c[i] << ' ';
	}
	cout << endl;
	Codeword<int> c2 = p.permutate(c);
	
	for (int i=0; i<c2.size(); i++) {
		cout << c2[i] << ' ';
	}
	cout << endl;
	Codeword<int> c3 = p.depermutate(c2);
	
	for (int i=0; i<c3.size(); i++) {
		cout << c3[i] << ' ';
	}
	cout << endl;
	return 0;
}