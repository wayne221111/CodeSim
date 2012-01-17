/*
 *  packet_receive_sim.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/5/19.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "randomc.h"
#include <ctime>
#define Run 100
using namespace std;

int main(){
	CRandomMersenne random(time(0));
	int K=10000, packetSize = 100;
	double inflation = 1.1, errorRate=0.06;
	unsigned long total=0;
	
	for (int run=0; run<Run; run++) {
		int received = 0, transmitted = 0;
		while (received < K*inflation) {
			transmitted += packetSize;
			if (random.Random() > errorRate) {
				received += packetSize;
			}
		}
		
		total += transmitted;
	}
	
	cout << total / ((double)(K)*Run);
	
	
}
