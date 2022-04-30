#ifndef GRID_FUNCTION_H
#define GRID_FUNCTION_H

#include "vectoroperators.hpp"
#include "../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>

namespace Function{
template <class T>
class UniformFunction1{
	double a,b;
	size_t N;
	double h;
	std::vector<T> F;
public:

    friend std::ostream& operator <<(std::ostream& os, const UniformFunction1<T>& F){
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.F.size();i++) {
            os << F.a + F.h*i << '\t' <<  debugdefs::to_debug_string(F.F[i])<< std::endl;
		}
		os << std::setprecision(prec);
		return os;
	}


	UniformFunction1(double _a = 0,double _b = 0,size_t _N = 0):a(_a),b(_b),N(_N),F(N){
		if(a>b)
			std::swap(a,b);
		h = (b-a)/(N-1);
	}
	UniformFunction1(size_t _N = 0):a(0),b(1),N(_N),F(_N){
		h = (b-a)/(N-1);
	}
	UniformFunction1(double a,double b,std::vector<T> Y):a(a),b(b),N(Y.size()),F(Y){
		if(a>b){
			std::swap(a,b);
			std::reverse(F.begin(),F.end());
		}
		h = (b-a)/(N-1);
	}
	
	template<typename FuncType>
	UniformFunction1(FuncType f,double _a = 0,double _b = 0,size_t _N = 0):a(_a),b(_b),N(_N),F(N){
		if(a>b)
			std::swap(a,b);
		h = (b-a)/(N-1);
		for(size_t i=0;i<N;++i){
			F[i] = f(a+h*i);
		}
	}
	std::vector<double> Xgrid() const{
		std::vector<double> XG(N);
		for(size_t i=0;i<N;i++){
			XG[i] = a + h*i;
		}
		return XG;
	}
	const std::vector<double> & Fval() const{
		return F;
	}
		
	template<typename FuncType>
	void map(FuncType f){
		for(size_t i=0;i<N;++i)
			F[i] = f(a+i*h);
	}
	
	size_t size()const {return N;}
	
	T operator ()(double x) const{
		int i = (x-a)/h;
		if(x < a)
			return F[0];
		else if(i >= N-1)
			return F[N-1];
		else{
			 double di = (x-a-i*h)/h;
			 return di*F[i+1]+(1-di)*F[i];
		}
	}
};
/*
template <typename T>
class UniformFunctionCubic1 : public UniformFunction1<T>{
	public:
	UniformFunctionCubic1();
	UniformFunctionCubic1(const UniformFunction1 &F):UniformFunction1(F){}
	T operator ()(double x) const{
		int i = (x-a)/h;
		if(x < a)
			return F[0];
		else if(i >= N-1)
			return F[N-1];
		else{
			if(i == 0)
				i = 1;
			else if(i >= N-2)
				i = N-3
				
			double x0 = (x-a-i*h)/h;
			double xm = xi+1.0;
			double xp = 1.0-xi;
			double xpp = 2.0-xi;
			
			return 0.5*xm*xpp*(F[i]*xp+F[i+1]*x0)-x0*xp/6*(F[i-1]*xpp+F[i+2]*xm);
		}
	}
};
*/
extern inline auto find_less(const std::vector<double> &X,double x){
	size_t N = X.size();
	size_t i1 = 0;
	size_t i2 = N-1;
	while(i1 +1 < i2){
		size_t i = (i1 + i2)/2;
		if(x < X[i]){
			i2 = i;
		}
		else {
			i1 = i;
		}
	}
	return i1;
}

template <typename T>
UniformFunction1<double> LoadUniformFunction1(const char *filename,bool isTitle = true){
	std::ifstream ifs(filename, std::ifstream::in);
	std::vector<double> Y;
	double xmin,xmax,x,y;
	if(isTitle){
		std::string S;
		std::getline(ifs,S);
	}
	bool first_time = true;
	while(ifs>>x>>y){
		if(first_time)
			xmin = x;
		xmax = x;
		Y.push_back(y);
	}
    return UniformFunction1<double>(xmin,xmax,Y);
}

template <class T>
class GridFunction1{
	std::vector<double> X;
	std::vector<T> Y;
	
	friend std::ostream& operator<<(std::ostream& os, const GridFunction1<T>& F){
		
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.size();i++) {
            os << F.X[i] << '\t' <<  debugdefs::to_debug_string(F.Y[i])<< std::endl;
		}
		os << std::setprecision(prec);
		return os;
	}
	
public:

	GridFunction1(size_t N = 0):X(N),Y(N){}
	enum SortState {SORTED,FLIPPED,CHAOTIC};
	
private:
	static void PrepareVectors(std::vector<double> &A,std::vector<T> &B,SortState SS){
		if(SS == SORTED)
			return;
		else if(SS == FLIPPED){
			std::reverse(A.begin(), A.end());
			std::reverse(B.begin(), B.end());
		}
		else{
			std::vector<std::pair<double,T>> AB(A.size());
			for(size_t i=0;i<AB.size();i++){
				AB[i] = std::pair<double,T>(A[i],B[i]);
			}
			std::sort(AB.begin(),AB.end(),
				[](const std::pair<double,T> &a,const std::pair<double,T> &b){
					return a.first<b.first;
				}
			);
			for(int i=0;i<AB.size();i++){
				A[i] = AB[i].first;
				B[i] = AB[i].second;
			}
		}
	}
	static void PrepareVector(std::vector<double> &A,SortState SS){
		if(SS == SORTED)
			return;
		else if(SS == FLIPPED)
			std::reverse(A.begin(), A.end());
		else{
			std::sort(A.begin(),A.end());
		}
	}
public:
	GridFunction1(const std::vector<double> &X):X(X),Y(X.size()){
		std::sort(this->X.begin(),this->X.end());
	}
	
	template <typename FuncType>
	GridFunction1(const std::vector<double> &X,FuncType f,bool sorted = false):
	X(X){
		PrepareVector(this->X);
		Y = vmap(f,this->X);
	}
	GridFunction1(const std::vector<double> &X,const std::vector<T> &Y,SortState SS = SORTED):
	X(X),Y(Y){
		PrepareVectors(this->X,this->Y,SS);
	}
	
	double &getX(size_t i){return X[i];}
	double &getY(size_t i){return Y[i];}
	
	const std::vector<double> &getXgrid() const{
		return X;
	}
	
	const std::vector<T> &getYgrid() const{
		return Y;
	}
	
	const size_t size()const{return X.size();}
	T operator ()(double x) const{
		size_t N = X.size();
		size_t i = find_less(X,x);
		size_t i1 = i+1;
		if(i1 >= N){
			return Y[i];
		}
		double a = (x-X[i])/(X[i1]-X[i]);
		return (1-a)*Y[i]+a*Y[i1];
	}
};

template <typename T>
std::vector<T> parse_string(const std::string &S){
	std::stringstream ss(S);
	std::vector<T> parsed;
	T it;
	while(!ss.eof()){
		ss >> it;
		parsed.push_back(it);
	}
	return parsed;
}

extern inline std::map<std::string, std::vector<double>> CSVTable(const std::string & filename){
	
	std::map<std::string, std::vector<double>> Funcs;
	
	std::ifstream ifs(filename, std::ifstream::in);
	
	if(!ifs.is_open()){
		throw std::runtime_error(std::string("CSVTable error: no such file: ") + filename);
	}
	
	std::string S;
	std::getline(ifs,S);
	std::vector<std::string> cols = parse_string<std::string>(S);
	
	for(auto title : cols){
		Funcs[title] = std::vector<double>();
	}
	
	std::vector<double> nums;
	while(!ifs.eof()){
		std::getline(ifs,S);
		nums = parse_string<double>(S);
		if(nums.size()>=cols.size()){
			for(size_t i=0;i<cols.size();i++){
				Funcs[cols[i]].push_back(nums[i]);
			}
		}
	}
	return Funcs;
};

};
#endif
