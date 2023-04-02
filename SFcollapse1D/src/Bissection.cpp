#include <iostream>
#include <cmath>

double bissection(double a, double b, double precision){
	double fine_tuning = 1.0;
	double pivot = std::abs((a-b)/b);
	int count = 0;
	while(fine_tuning > precision && count < 10000){
		count = count+1;
		fine_tuning =  pivot - cos(pivot);
		if(cos(a)*cos(pivot)<0){
			pivot =std::abs((a+pivot)/2);
		}else if(cos(b)*cos(pivot>0)){
			pivot = std::abs((b+pivot)/2);
		}else{break;
			std::cout << "pivot = raiz"<< std::endl;}
	}
	return pivot;
}

int main(){
	//std::cout << bissection(1,2,0.00001) << std::endl;
	
return 0;
}

/*
vou testar inicialmente com a função x = cos(x)
*/



