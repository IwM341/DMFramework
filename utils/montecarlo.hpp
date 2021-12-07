#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <random>

namespace MC{
	
template <typename T>
struct MCResult{
	T Result;
	double RemainDensity;
	MCResult(T Result,double RemainDensity):Result(Result),RemainDensity(RemainDensity){}
};

};
#endif
