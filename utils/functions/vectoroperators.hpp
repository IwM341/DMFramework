#ifndef VECTOR_OPERATORS_H
#define VECTOR_OPERATORS_H

#include <vector>
#include <type_traits>
#include <utility>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
//PRINT VECTOR

struct empty_type{
    template <typename T>
    empty_type & operator = (const T &x){return *this;}
};
static empty_type _;

#define _T std::make_tuple
#define _R make_ref_typle

template <typename...Args>
auto make_ref_typle(Args&...args){
    return std::tuple<Args&...>(args...);
}

template <typename V>
inline auto iterator_last(V & list){
    return list.begin() + (list.size()-1);
}
template <typename V>
inline auto const_iterator_last(const V & list){
    return list.cbegin() + (list.size()-1);
}

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

template <typename Container>
bool contains_nan(const Container&C){
    for(const auto & val : C){
        if(std::isnan(val)){
            return true;
        }
    }
    return false;
}

template <class T>
auto vector_reverse(const std::vector<T> &X){
    std::vector<T> Y = X;
    std::reverse(Y.begin(),Y.end());
    return Y;
}
template <class T>
auto vector_reverse(std::vector<T> &&X){
    std::reverse(X.begin(),X.end());
    return X;
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


template <typename T>
T max(const std::vector<T> &X){
    T tmp = X[0];
    for(size_t i=1;i<X.size();++i){
        if(tmp < X[i])
            tmp = X[i];
    }
    return tmp;
}
template <typename T>
T min(const std::vector<T> &X){
    T tmp = X[0];
    for(size_t i=1;i<X.size();++i){
        if(tmp > X[i])
            tmp = X[i];
    }
    return tmp;
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
    if(!V.size())
        return 0;
    T summ = V[0];
    for(size_t i=1;i<V.size();++i)
        summ += V[i];
    return summ;
}


template <typename Functype>
inline auto sum_lambda(size_t size,Functype && F){
    typename std::decay<decltype(F(0))>::type summ = 0;
    for(size_t i=0;i<size;++i)
        summ += F(i);
    return summ;
}

//CREATE VECTOR FROM LAMBDA
template<typename LambdaType>
inline auto Vector(size_t N, LambdaType &&F){
    std::vector<decltype(F(0))> Y(N);
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
    std::vector<std::remove_reference_t<decltype(std::declval<T1>()-std::declval<T2>())>> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]-Y[i];
	return Z;
}
//vector - value
template <typename T1,typename T2>
inline auto operator - (const std::vector<T1> & X,const T2 & Y){
    std::vector<std::remove_reference_t<decltype(std::declval<T1>()-std::declval<T2>())>> Z(X.size());
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
    std::vector<T> Z(Y.size());
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

template <typename IteratorType,typename T>
void __fill__iterator(IteratorType start,const T &first){
    *start = first;
    ++start;
}

template <typename IteratorType,typename T,typename...Args>
void __fill__iterator(IteratorType start,const T &first,const Args&...Other){
    *start = first;
    ++start;
    __fill__iterator(start,Other...);
}

template <typename...Args>
struct arg_count;

template <typename T,typename...Args>
struct arg_count<T,Args...>{
  constexpr static size_t N = 1 + arg_count<Args...>::N;
};
template <>
struct arg_count<>{
  constexpr static size_t N = 0;
};
#endif
