/*
 *  Byte.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/5/5.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */



#include<vector>
#include<string>
#include"CodeSim.h"
using namespace std;

namespace CodeSim {
	
	
	
	
	Byte::Byte(){
		value = 0;
		setErased(false);
	}
	
	//template<class T>
	Byte::Byte(int t){
		if (t == -1) {
			setErased(true);
		}
		else{
			value = t & mask;
			setErased(false);
		}
	}
	
	//	template<class T>
	//	Byte Byte::operator=(T t){
	//		value = t;
	//		return *this;
	//	}
	Byte Byte::operator+(Byte t){
		return Byte((value ^ t.value) );
	}
	
	//	Byte Byte::operator*(Byte t){
	//		return Byte(this->value & t.value);
	//	}
	
	string Byte::toString()
	{
		if (isErased()) {
			return "********";
		}
		unsigned int tmp = value;
		string s = "";
		for (int i=0; i<8; i++) {
			if (tmp & 1) {
				s = '1'+s;
			}
			else {
				s = '0'+s;
			}
			
			tmp = tmp>>1;

		}
		return s;
	}
	
	
}
