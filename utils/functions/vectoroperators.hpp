#ifndef VECTOR_OPERATORS_H
#define VECTOR_OPERATORS_H

#include <vector>
#include <type_traits>
#include <utility>
#include <iostream>
#include <sstream>
//PRINT VECTOR



namespace std{

template <class T>
std::ostream & operator << (std::ostream & os,const std::vector<T>& V){
    os << "Vector[ ";
    for (size_t i=0;i<V.size();i++) {
        if(i < V.size()-1)
            os << V[i] << ", ";
        else
            os << V[i] << "]";
    }
    return os;
}


template <class T>
std::string to_string(const std::vector<T>& V){
    std::ostringstream os;
    os << V;
    return os.str();
}
};

//MAP FUNCTION ON VECTOR
template <typename T,typename FuncType>
inline auto vmap(const FuncType &F, const std::vector<T> &X){
	std::vector<typename std::invoke_result<FuncType,T>::type> Y(X.size());
	for(size_t i=0;i<X.size();i++){
		Y[i] = F(X[i]);
	}
	return Y;
}

template <typename IteratorType>
inline size_t unique_values_sorted(IteratorType start,IteratorType end){
    size_t sum = 0;
    if(start == end)
        return sum;

    auto tmp_value = *start;
    sum++;
    for(;start!=end;++start){
        if(*start != tmp_value){
            sum++;
            tmp_value = *start;
        }
    }
    return sum;
}
/*
template <typename T>
inline size_t unique_values_unsorted(const std::vector<T> & V){
    size_t sum = 0;
    if(!V.size())
        return sum;

    T tmp_value = V[0];
    sum++;
    for(size_t i = 1;i<V.size();++i){
        if(V[i] != tmp_value){
            sum++;
            tmp_value = V[i];
        }
    }
    return sum;
}
*/
template <typename T>
inline T vector_sum(const std::vector<T> & V){
    T summ = 0;
    for(auto it : V)
        summ += it;
    return summ;
}
//CREATE VECTOR FROM LAMBDA
template<typename LambdaType>
inline auto Vector(size_t N,const LambdaType &F){
	std::vector<typename std::invoke_result<LambdaType, int>::type> Y(N);
	for(size_t i=0;i<N;++i){
		Y[i] = F(i);
	}
	return Y;
}


//SUM VECTORS
//vector + vector
template <typename T1,typename T2>
inline auto operator + (const std::vector<T1> & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()+std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]+Y[i];
	return Z;
}
//vector + value
template <typename T1,typename T2>
inline auto operator + (const std::vector<T1> & X,const T2 & Y){
	std::vector<decltype(std::declval<T1>()+std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]+Y;
	return Z;
}
//value + vector
template <typename T1,typename T2>
inline auto operator + (const T1 & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()+std::declval<T2>())> Z(Y.size());
	for(size_t i=0;i<Y.size();i++)
		Z[i] = Y[i]+X;
	return Z;
}

//SUBSTRACTION VECTORS
//vector - vector
template <typename T1,typename T2>
inline auto operator - (const std::vector<T1> & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()-std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]-Y[i];
	return Z;
}
//vector - value
template <typename T1,typename T2>
inline auto operator - (const std::vector<T1> & X,const T2 & Y){
	std::vector<decltype(std::declval<T1>()-std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]-Y;
	return Z;
}
//value - vector
template <typename T1,typename T2>
inline auto operator - (const T1 & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()-std::declval<T2>())> Z(Y.size());
	for(size_t i=0;i<Y.size();i++)
		Z[i] = X - Y[i];
	return Z;
}
//- vector
template <typename T>
inline auto operator -(const std::vector<T> & Y){
	std::vector<decltype(std::declval<T>())> Z(Y.size());
	for(size_t i=0;i<Y.size();i++)
		Z[i] = - Y[i];
	return Z;
}

//ADD ASSIGNMENT
//vector += vector
template <typename T1,typename T2>
inline const auto & operator += (std::vector<T1> & X,const std::vector<T2> & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]+=Y[i];
	return X;
}
//vector += value
template <typename T1,typename T2>
inline const auto & operator += (std::vector<T1> & X,const T2 & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]+=Y;
	return X;
}

//SUBSTRACTION ASSIGNMENT
//vector -= vector
template <typename T1,typename T2>
inline const auto & operator -= (std::vector<T1> & X,const std::vector<T2> & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]-=Y[i];
	return X;
}
//vector -= value
template <typename T1,typename T2>
inline const auto & operator -= (std::vector<T1> & X,const T2 & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]-=Y;
	return X;
}



//MULTIPLYING VECTORS
//vector * vector
template <typename T1,typename T2>
inline auto operator * (const std::vector<T1> & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()*std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]*Y[i];
	return Z;
}
//vector * value
template <typename T1,typename T2>
inline auto operator * (const std::vector<T1> & X,const T2 & Y){
	std::vector<decltype(std::declval<T1>()*std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]*Y;
	return Z;
}
//value * vector
template <typename T1,typename T2>
inline auto operator * (const T1 & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()*std::declval<T2>())> Z(Y.size());
	for(size_t i=0;i<Y.size();i++)
		Z[i] = Y[i]*X;
	return Z;
}

//DIVISION VECTORS
//vector / vector
template <typename T1,typename T2>
inline auto operator / (const std::vector<T1> & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()/std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]/Y[i];
	return Z;
}
//vector - value
template <typename T1,typename T2>
inline auto operator / (const std::vector<T1> & X,const T2 & Y){
	std::vector<decltype(std::declval<T1>()/std::declval<T2>())> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]/Y;
	return Z;
}
//value - vector
template <typename T1,typename T2>
inline auto operator / (const T1 & X,const std::vector<T2> & Y){
	std::vector<decltype(std::declval<T1>()/std::declval<T2>())> Z(Y.size());
	for(size_t i=0;i<Y.size();i++)
		Z[i] = X / Y[i];
	return Z;
}


//MULT ASSIGNMENT
//vector *= vector
template <typename T1,typename T2>
inline const auto & operator *= (std::vector<T1> & X,const std::vector<T2> & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]*=Y[i];
	return X;
}
//vector *= value
template <typename T1,typename T2>
inline const auto & operator *= (std::vector<T1> & X,const T2 & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]*=Y;
	return X;
}

//DIVISION ASSIGNMENT
//vector /= vector
template <typename T1,typename T2>
inline const auto & operator /= (std::vector<T1> & X,const std::vector<T2> & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]/=Y[i];
	return X;
}
//vector -= value
template <typename T1,typename T2>
inline const auto & operator /= (std::vector<T1> & X,const T2 & Y){
	for(size_t i=0;i<X.size();i++)
		X[i]/=Y;
	return X;
}
#endif
