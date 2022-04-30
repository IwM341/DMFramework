#ifndef COMPLEX_EX_H
#define COMPLEX_EX_H
#include <complex>
const std::complex<double> Id(0,1);

template<typename U>
auto operator +(std::complex<U> x,int y){
    return std::complex<U>(x.real()+y,x.imag());
}
template<typename U>
auto operator +(int y,std::complex<U> x){
    return std::complex<U>(x.real()+y,x.imag());
}

template<typename U>
auto operator -(std::complex<U> x,int y){
    return std::complex<U>(x.real()-y,x.imag());
}
template<typename U>
auto operator -(int y,std::complex<U> x){
    return std::complex<U>(y-x.real(),x.imag());
}


template<typename U>
auto operator *(std::complex<U> x,int y){
    return std::complex<U>(x.real()*y,x.imag()*y);
}
template<typename U>
auto operator *(int y,std::complex<U> x){
    return std::complex<U>(x.real()*y,x.imag()*y);
}

template<typename U>
auto operator /(std::complex<U> x,int y){
    return std::complex<U>(x.real()/y,x.imag()/y);
}
template<typename U>
auto operator /(int y,std::complex<U> x){
    auto c = y/(x.real()*x.real()+x.imag()*x.imag());
    return std::complex<U>(x.real()*c,-x.imag()*c);
}

template<typename U>
U fast_pow(U x,int n){
    if(n<0){
        x = static_cast<U>(1)/x;
        n = -n;
    }
    U res = static_cast<U>(1);
    while(n > 0){
        if(n%2)
            res *= x;
        x = x*x;
        n = n/2;
    }
    return res;
}


#endif
