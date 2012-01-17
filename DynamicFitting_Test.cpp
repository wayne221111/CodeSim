/*
 *  DynamicFitting_Test.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/12/16.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "DynamicFitting.h"
#include <vector>
using namespace std;

int main(){
	const int NPT=26,MA=4;
	const DP gues_d[MA]={0,0,0,0};
	DP data_x[NPT] = {0.0000,0.0100,0.0200,0.0300,0.0400,0.0500,
		0.0600,0.0700,0.0800,0.0900,0.1000,0.1100,0.1200,0.1300,
		0.1400,0.1500,0.1600,0.1700,0.1800,0.1900,0.2000,0.2100,
		0.2200,0.2300,0.2400,0.2500}
	, data_y[NPT] = {-0.2132,-0.2236,-0.2590,-0.2872,-0.3388,
		-0.3951,-0.4836,-0.5917,-0.7452,-0.9313,-1.1605,
		-1.4161,-1.7149,-2.0230,-2.3408,-2.6625,-3.0022,
		-3.2591,-3.4218,-3.5401,-3.6563,-3.7474,-3.8143,
		-3.8747,-3.9219,-3.9629};
	vector<double> x(data_x, data_x+NPT), y(data_y, data_y+NPT), a(gues_d, gues_d+MA);
	double ss;
	
	//while (cin >> a[0]) 
	{
		//for (int i = 1; i<MA; i++) {
//			cin >> a[i];
//		}
		int k = regression(x, y, a, cubic,ss);
		
		cout << ss;
		for (int i=0; i<MA; i++) {
			cout << '\t' << a[i];
		}
		
		cout << '\n';
	}
	
	//dynamicFitting(x,y,sigmoid,sigmoidRange,a,ss);
//	cout << "Final:\nargs:";
//	for (int i=0; i<MA; i++) {
//		cout << ' ' << a[i];
//	}
//	cout << "\nSS: " << ss;
	//cout << "\nAfter " << k << " iterations.\n";
}