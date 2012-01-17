/*
 *  Symbol.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/27.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include<vector>
#include<string>
#include"CodeSim.h"
using namespace std;

namespace CodeSim {

	
	
	
	Bit::Bit(){
		value = 0;
		setErased(false);
	}
		
	//template<class T>
	Bit::Bit(int t){
		if (t == -1) {
			setErased(true);
		}
		else{
			value = t;
			setErased(false);
		}
	}
		
//	template<class T>
//	Bit Bit::operator=(T t){
//		value = t;
//		return *this;
//	}
	Bit Bit::operator+(Bit t){
		return Bit(value ^ t.value);
	}
	
//	Bit Bit::operator*(Bit t){
//		return Bit(this->value & t.value);
//	}
	
	string Bit::toString()
	{
		if (isErased()) {
			return "*";
		}
		if (value) {
			return "1";
		}
		else {
			return "0";
		}
	}
		
	
}
