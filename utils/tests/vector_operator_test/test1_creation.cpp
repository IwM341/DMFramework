#include <vectoroperators.hpp>
#include <iostream>
int main(void){
	auto V = Vector(10,[](size_t i){return i*i+0.1;});
	std::cout << V <<std::endl;
	std::cout << vmap([](double i){return 0.5*(i+1.0/i);},V) <<std::endl;
	return 0;
}
