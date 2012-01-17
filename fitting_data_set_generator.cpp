/*
 *  fitting_data_set_generator.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/11/14.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;



int main(){
	string s;
	vector< vector<double> > data;

	
	while(getline(cin, s)){
		data.clear();
		{
			istringstream iss(s);
			double d;
			int i=0;
			while (iss >> d) {
				data.push_back(vector<double>());
				data[i].push_back(d);
				i++;
			}
			if (i==0) {
				cout << '\n';
				continue;
			}
		}
		
		while ( getline(cin, s) ) {
			istringstream iss(s);
			double d;
			int i;
			for (i=0; iss >> d; i++) {
				data[i].push_back(d);
			}
			if (i==0) {
				break;
			}
		}
		
		for (vector< vector<double> >::iterator i=data.begin(); i!=data.end(); i++) {
			
			sort(i->begin(), i->end());
			*i = vector<double>(i->rbegin(), i->rend());
			
			int size = i->size();
			for (int s=0; s<100; s++) {
				double sum=0;
				int n=0;
				for (int x=s*(size/100); x< (s+1)*(size/100); x++) {
					sum += (*i)[x];
					n++;
				}
//				if (sum > 0) 
					cout << sum/n << '\t';
//				else {
//					cout << "1.00E-8\t";
//				}

			}
			cout << '\n';
		}
		
		cout << '\n';
		
	}
	return 0;
}