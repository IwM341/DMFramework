#ifndef INT_LAMBDA
#define INT_LAMBDA

//#include "omp.h"
#include <functional>
#include <utility>
#include <cmath>
#include <vector>
	
template <typename FuncType>			
extern inline double integrateAB(FuncType f,double a,double b,int N){
	
	double sum = 0.0;
	double h = (b-a)/N;
	for (int i = 0; i < N; i++) {
		sum+= f(a + (i+0.5)*h)*h;
	}
	return sum;
}

template <typename FuncType>			
extern inline double integrateABP(FuncType f,double a,double b,int N){
	
	double sum = 0.0;
	double h = (b-a)/N;
	
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < N; i++) {
		sum += f(a + (i+0.5)*h)*h;
	}
	return sum;
}

const double Xh[5] = {0.5, 0.5+0.5*.53846931010568309105, 0.5-0.5*.53846931010568309105,
					0.5+0.5*.90617984593866399282, 0.5-0.5*.90617984593866399282};
const double Wh[5] = {.28444444444444444444, .23931433524968323402, .23931433524968323402,
					.11846344252809454376, .11846344252809454376};

#define X0 0.5
#define X1 0.76923465505284154552
#define X2 0.23076534494715845448
#define X3 0.95308992296933199641
#define X4 0.04691007703066800359

#define W0 0.28444444444444444444
#define W1 0.23931433524968323402
#define W2 0.23931433524968323402
#define W3 0.11846344252809454376
#define W4 0.11846344252809454376


template <typename FuncType>	
extern inline double integrateAB5(FuncType f,double a,double b,int N){

	double sum = 0.0;
	double xi;
	double h = (b-a)/N;
	double h0 = h*X0,h1 = h*X1,h2 = h*X2, h3 = h*X3, h4 = h*h4;
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += h*(W0*f(xi + h0)+
					W1*f(xi + h1)+
					W2*f(xi + h2)+
					W3*f(xi + h3)+
					W4*f(xi + h4));//sum5*h;
	}
	return sum;
}

template <typename FuncType>	
extern inline double integrateAB5P(FuncType f,double a,double b,int N){

	double sum = 0.0;
	double xi;
	double h = (b-a)/N;
	double h0 = h*X0,h1 = h*X1,h2 = h*X2, h3 = h*X3, h4 = h*h4;
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += h*(W0*f(xi + h0)+
					W1*f(xi + h1)+
					W2*f(xi + h2)+
					W3*f(xi + h3)+
					W4*f(xi + h4));//sum5*h;
	}
	return sum;
}

template <typename FuncType>
extern inline double integrateLeftLogAB5(FuncType f,double a,double b,double eps,const int N){
	const double L = b-a;
	
	const double T = pow(eps/(L+eps),1.0/N);
	
	double sum = 0.0;
	double sum5;
	
	

	//#pragma omp parallel for reduction(+:sum) shared(N,T,L,eps)
	for (int i = 0; i < N; i++) {
		//sum5 = 0.0;
		const double xi = a - eps + (L+eps)*pow(T,i+1);
		const double xi1 = a - eps + (L+eps)*pow(T,i);
		const double h = xi1 - xi;
		//std::cout<<xi <<" , " << xi1 << " , " << h <<std::endl;
		//for(int j = 0;j<5;j++)
			//sum5 += Wh[j]*f(xi + h*Xh[j]);
		sum += h*(W0*f(xi + h*X0)+
					W1*f(xi + h*X1)+
					W2*f(xi + h*X2)+
					W3*f(xi + h*X3)+
					W4*f(xi + h*X4));//sum5*h;
	}
	return sum;
}	

template <typename FuncType>
extern inline double integrateRightLogAB5(FuncType f,double a,double b,double eps,const int N){
	const double L = b-a;
	
	const double T = pow(eps/(L+eps),1.0/N);
	
	double sum = 0.0;
	double sum5;
	
	
	//#pragma omp parallel for reduction(+:sum)  shared(N,T,L,eps)
	for (int i = 0; i < N; i++) {
		//sum5 = 0.0;
		const double xi = b + eps - (L+eps)*pow(T,i);
		const double xi1 = b + eps - (L+eps)*pow(T,i+1);
		const double h = xi1 - xi;
		//for(int j = 0;j<5;j++)
			//sum5 += Wh[j]*f(xi + h*Xh[j]);
		sum += h*(W0*f(xi + h*X0)+
					W1*f(xi + h*X1)+
					W2*f(xi + h*X2)+
					W3*f(xi + h*X3)+
					W4*f(xi + h*X4));//sum5*h;
	}
	return sum;
}	

template <typename T>
extern inline double integrateAB2(std::vector<double> X,std::vector<T> F){
	T sum = 0;
	int N = X.size();
	for (int i = 0; i < N-1; i++){
		sum += 0.5*(F[i]+F[i+1])*(X[i+1]-X[i]);
	}
	return sum;
}

			
#endif
