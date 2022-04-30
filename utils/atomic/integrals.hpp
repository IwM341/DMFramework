#ifndef INTEGRALS_HPP
#define INTEGRALS_HPP
#include <complex>
#include <cmath>
#include "../complex/complex_ex.hpp"
#include "../functions/vectoroperators.hpp"

namespace intalg {

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

template <typename U>
std::complex<U> PI0k(int k,U v){
    if(k == 1){
        return 2*atan(v)/v;
    }
    else{
        std::complex<U> z (0,v);
        return(1.0/(z*(k-1)*fast_pow(1.0-z,k-1)) - 1.0/(z*(k-1)*fast_pow(1.0+z,k-1)));
    }
}
template <typename U>
std::complex<U> PI1k(int k, U v){
    if(k == 0){
        return 0;
    }
    return (PI0k(k-1,v)-PI0k(k,v))/(Id*v);
}

template <typename U>
std::complex<U> PIl0(int l,U v){
    return (l == 1? 2.0 : 0.0);
}



template <typename U1,typename U2,typename U3,typename U>
std::complex<U> PI_Next(int l,/*P(l-2,k) = */U1 PI_l_2_k,
                           /*P(l-1,k-1) = */U2 PI_l_1_k_1,
                           /*P(l-1,k) = */U3 PI_l_1_k, U v)
{
    return std::complex<double>(0,(2*l-1)/(v*l))*
                (PI_l_1_k-PI_l_1_k_1)-
            (l-1.0)/l*PI_l_2_k;
}

template <typename U>
std::complex<U> PI_Next2(int l_tmp, int k, std::complex<U> PI_l_tmp_k_2,
                           std::complex<U> PI_l_tmp_k_1, U v)
{
    if(k<3){
        std::cout << "error in PI_Next2: k<3, division on 0" << std::endl;
        return 1.0/0.0;
    }
    return (2*(k-2)*(k-2)*PI_l_tmp_k_1 + (l_tmp*(l_tmp+1) - (k-2)*(k-3))*PI_l_tmp_k_2)/((k-1)*(k-2)*(1+v*v));
}

template <typename U>
/*the k from 0 to MaxK*/
std::vector<std::complex<U>> PIl0l1Vec(bool l,int MaxK, U v){
    if(l){
        return Vector(MaxK+1,[v](size_t i){
            return PI1k(i,v);
        });
    }
    else{
        return Vector(MaxK+1,[v](size_t i){
            return PI0k(i,v);
        });
    }
};

template <typename U>
void PI_Next(const std::vector<std::complex<U>> &PI_l_2,
             const std::vector<std::complex<U>> &PI_l_1,
             std::vector<std::complex<U>> &Result,int l, U v){
    if(!(Result.size() == PI_l_2.size() && PI_l_2.size() == PI_l_1.size())){
        std::cout << "error at PI_Next(): PI_l_2.size() = "<<PI_l_2.size()<<
                     ", PI_l_1.size() = "<<PI_l_1.size()<<
                     ", Result.size() = "<<Result.size()<<" -  should be the same" <<std::endl;
        return;
    }
    if(!Result.size()){
        return;
    }
    else{
        Result[0] = 0;
    }
    for(int i=1;i<Result.size();++i){
        Result[i] = PI_Next(l,PI_l_2[i],PI_l_1[i-1],PI_l_1[i],v);
    }
}


template <typename U>
auto PI_Next(const std::vector<std::complex<U>> &PI_l_2,
             const std::vector<std::complex<U>> &PI_l_1,int l_next, U v){
    std::vector<std::complex<U>> PI_Res(PI_l_2.size());
    PI_Next(PI_l_2,PI_l_1,PI_Res,l_next,v);
    return PI_Res;
}

template <typename U>
auto PI_Table(int Lmax,int Kmax,U v){
    std::vector<std::vector<std::complex<U>>> PI_Tab(Lmax+1);
    if(Lmax >= 0){
        PI_Tab[0] = PIl0l1Vec(0,Kmax,v);
    }
    if(Lmax > 0){
        PI_Tab[1] = PIl0l1Vec(1,Kmax,v);
    }
    for(int l=2;l<=Lmax;++l){
        PI_Tab[l] = PI_Next(PI_Tab[l-2],PI_Tab[l-1],l,v);
    }
    return PI_Tab;
}

template <typename U>
auto PI_l_Vec(int l,int Kmax, U v){
    std::vector<std::complex<U>> Result(Kmax+1);
    if(Kmax >= 0)
        Result[0] = PIl0(l,v);
    if(Kmax > 0){
        std::complex<U> Rlm1 = PI0k(1,v);
        std::complex<U> Rlm2 = PI0k(2,v);
        std::complex<U> Rlp1 = PI1k(1,v);
        std::complex<U> Rlp2 = PI1k(2,v);

        for(int ltmp = 2;ltmp<=l;++ltmp){
            auto Rl1tmp = PI_Next(ltmp,Rlm1,0.0,Rlp1,v);
            auto Rl2tmp = PI_Next(ltmp,Rlm2,Rlp1,Rlp2,v);
            Rlm1 = Rlp1;
            Rlm2 = Rlp2;
            Rlp1 = Rl1tmp;
            Rlp2 = Rl2tmp;
        }
        Result[1] = Rlp1;
        if(Kmax > 1){
            Result[2] = Rlp2;
        }
        for(int i=3;i<=Kmax;++i){
            Result[i] = PI_Next2(l,i,Result[i-2],Result[i-1],v);
        }
    }
    return Result;
}

};

#endif // INTEGRALS_HPP
