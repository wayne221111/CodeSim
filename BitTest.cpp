/*
 *  BitTest.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/26.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include<iostream>
#include"CodeSim.h"

using namespace std;
using namespace CodeSim;

int main(){
	Bit a = 3, b = 0;
	cout << "Size of Bit: " << sizeof(Bit) << endl;
	cout << "a+b = " << (a+b).toString() << endl;
	cout << "a*b = " << (a*b).toString() << endl;
	cout << "a+a = " << (a+a).toString() << endl;
	cout << "a*a = " << (a*a).toString() << endl;
	cout << "b+b = " << (b+b).toString() << endl;
	cout << "a+1 = " << (a+1).toString() << endl;
}
