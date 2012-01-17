/*
 *  histogram_cumulator.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/6/2.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <map>

using namespace std;

int main(){
	map<long, unsigned long> histo;
	long a;
	unsigned long b, sum=0;
	
	while (cin >> a) {
		if(cin >> b){
			histo[a] += b;
			sum += b;
		}
	}
	
	cout << "\n\n\nkey\tcount\n";
	for (map<long, unsigned long>::iterator i=histo.begin(); i!=histo.end(); i++) {
		cout << i->first << '\t' << i->second << '\n';
	}
	cout << "total:\t" << sum << "\n\n";
	
//	cout << "key\tfrequency\n";
//	for (map<long, unsigned long>::iterator i=histo.begin(); i!=histo.end(); i++) {
//		cout << i->first << '\t' << i->second/(double) sum << '\n';
//	}
	
	return 0;
}