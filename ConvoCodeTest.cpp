/*
 *  ConvoCode.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "ConvoCode.h"
#include "LTCode.h"
#include "randomc.h"
#include<time.h>
#define L 10000000

using namespace std;
using namespace CodeSim;
int main(){
	CRandomMersenne r(time(0));
	int startTime = time(0);
	cout << octToDec(19) << "\n";
	ConvoCode cc("convo.txt");
	cc.showInfo();
	
	cout << "K: " << cc.getK()<< " N: " << cc.getN() << endl;
	Codeword<Bit> a;
	for (int i=0; i<2*L; i++) {
		a.push_back(r.IRandomX(0, 1));
	}
	
//	cout << "input:  ";
//	for (int i = 0; i<a.size(); i++) {
//		cout << a[i].toString();
//	}
//	cout << '\n';
	cout << "encoding...\n";
	Codeword<Bit> b = cc.encode(a);
	
	for (int i = 0; i<b.size(); i++) {
		if (r.Random() < 0.1) {
			b[i].setErased(true);
		}
		//cout << b[i].toString();
	}
	
	
	cout << "decoding...\n";
	Codeword<Bit> b2 = cc.decode(b);
	int sum[2] = {0,0};
	
	
	for (int i = 0; i<b2.size(); i++) {
		if (b2[i].getValue() != a[i].getValue()) {
			sum[i%2]++;
		}
	}
	cout << '\n';
	
	cout << sum[0] / (double)L << ' ' << sum[1] / (double)L << '\n';
	cout << "Time: " << time(0) - startTime << '\n';
	return 0;
}