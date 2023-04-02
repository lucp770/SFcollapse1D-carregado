// 
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

double funcao(double x){
	double result = x*x - 2*x ;
	return result;
}

bool is_inclination_positive(double x){

	double delta = 0.0000001;
	double x1 = x+delta;
	double y0 = x*x - 2*x;
	double y1 = x1*x1 - 2*x1;
	double inclination = (y1-y0)/(x1-x);

	if(inclination > 0){
		return true;
	}else{
		return false;
	}
}


int main(){

	double x = 2.0;
	cout << to_string(is_inclination_positive(x)) << endl;

	return 0;
}
