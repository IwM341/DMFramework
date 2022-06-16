#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <random>
#include <utility>
namespace MC{
	
template <typename T>
struct MCResult{
	T Result;
	double RemainDensity;
	MCResult(T Result,double RemainDensity):Result(Result),RemainDensity(RemainDensity){}
};


template <typename T>
struct MCIntegral{
    T Result;
    T Sigma;
    MCIntegral(T Result,T Sigma = 0):Result(Result),Sigma(Sigma){}

    operator T (){
        return Result;
    }
    operator std::pair<T,T> (){
        return std::pair<T,T>(Result,Sigma);
    }
};

template<typename GenType>
auto MCIntegrate(GenType F,int N){
    double sum = 0;
    double sum2 = 0;
    for (size_t i=0;i<N;++i){
        auto f = F();
        sum += f;
        sum2 += f*f;
    }
    sum = sum/N;
    sum2 = sum2/N;
    return MCIntegral<double>(sum,sqrt(sum2-sum*sum));
}


};
#endif
