/*
 *  NR_test.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/12/12.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <vector>
#include "nr.h"
using namespace std;


void sigmoid(const DP x, Vec_I_DP &a, DP &y, Vec_O_DP &dyda)
{
	int i;
	DP fac,ex,arg;
	
	if(a.size()!=5)
		cerr << "sigmoid: wrong # of parameters\n";
	ex = pow(1+exp((-x+a[3])/a[1]),(-a[2]));
	y=a[4]+a[0]*ex;
	dyda[0] = ex;
	dyda[2] = -a[0]*log(1+exp((-x+a[3])/a[1]))*ex;
	ex /= 1+exp((-x+a[3])/a[1]);
	dyda[1] = a[0]*a[2]*(-x+a[3])*exp((-x+a[3])/a[1])/a[1]/a[1]*ex;
	dyda[3] = -a[0]*a[2]*exp((-x+a[3])/a[1])/a[1]*ex;
	dyda[4] = 1;
	//	dyda[0] = 1/(pow(1+exp(-(x-a[3])/a[1]),a[2]));
	//	dyda[1] = -a[0]/(pow(1+exp(-(x-a[3])/a[1]),a[2]))*a[2]*(x-a[3])/(a[1]*a[1])
	//		*exp(-(x-a[3])/a[1])/(1+exp(-(x-a[3])/a[1]));
	//	dyda[2] = -a[0]/(pow(1+exp(-(x-a[3])/a[1]),a[2]))*log(1+exp(-(x-a[3])/a[1]));
	//	dyda[3] = -a[0]/(pow(1+exp(-(x-a[3])/a[1]),a[2]))*a[2]/a[1]*exp(-(x-a[3])/a[1])/(1+exp(-(x-a[3])/a[1]));
}

int main(){
	vector<DP> x;
	double t;
	while ( cin >> t) {
		if (t<0) {
			break;
		}
		x.push_back(t);
		cout << t<<endl;
	}
	DP arg[] = {3.8179, -0.0301, 1.1329, 0.1362, -3.9848};
	Vec_I_DP a(arg, 5);
	double y;
	Vec_O_DP dyda(0.0,5);
	for (int i=0; i<x.size(); i++) {
		sigmoid(x[i], a, y, dyda);
		cout << y << endl;
	}
}
