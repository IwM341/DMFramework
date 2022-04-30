#ifndef POLYNOMEXTEND_H
#define POLYNOMEXTEND_H

#include <boost/math/tools/polynomial.hpp>
#include <math.h>
namespace mtls = boost::math::tools;


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
}

template <typename U,typename V>
U IntegratePolynom(const mtls::polynomial<U> & P, V a,V b){
    auto P1 = P.integrate();
    return P1(b)-P1(a);
}

template <typename U>
class extended_polynomial : public mtls::polynomial<U>{
    int deg;//nagative degree
public:
    typedef typename std::vector<U>::value_type value_type;
    typedef typename std::vector<U>::size_type size_type;

    template <typename V>
    extended_polynomial (const mtls::polynomial<V> &P,int deg = 0):mtls::polynomial<U>(P),deg(deg){}

    extended_polynomial (const mtls::polynomial<U> &&P, int deg =0):mtls::polynomial<U>(std::move(P)),deg(deg){}

    template <typename V>
    extended_polynomial(const V* data, unsigned order, int deg =0):mtls::polynomial<U>(data,order),deg(deg){}

    template <class I>
    extended_polynomial(I first, I last, int deg = 0):mtls::polynomial<U>(first,last),deg(deg){}

    extended_polynomial(std::vector<U>&& p,int deg =0):mtls::polynomial<U>(std::move(p)),deg(deg){}
    //antipolynomial():mtls::polynomial<U>(){}

    template <class Range, typename std::enable_if<boost::math::tools::detail::is_const_iterable<Range>::value, bool>::type = true>
    explicit extended_polynomial(const Range& r,int deg =0):mtls::polynomial<U>(r),deg(deg){}

    extended_polynomial(std::initializer_list<U> l,int deg =0):mtls::polynomial<U>(l),deg(deg){}

    extended_polynomial(int up_deg,int down_deg = 0):
        mtls::polynomial<U>(std::vector<U>(up_deg-down_deg+1,0)),deg(down_deg){}

    extended_polynomial & operator =(extended_polynomial &&p){
        mtls::polynomial<U>::operator=(p);
        return *this;
    }
    extended_polynomial & operator =(const extended_polynomial &p){
        mtls::polynomial<U>::operator=(p);
        return *this;
    }

    U operator ()(U z) const {
        return mtls::polynomial<U>::operator()(z)*fast_pow(z,deg);
    }
    inline size_type up_degree() const{
        return this->degree()+deg;
    }
    inline size_type down_degree() const{
        return deg;
    }

    value_type& operator[](size_type i)
    {
       return mtls::polynomial<U>::m_data[i-deg];
    }

    const value_type& operator[](size_type i) const
    {
       return mtls::polynomial<U>::m_data[i-deg];
    }

    operator mtls::polynomial<U> () const{
        std::vector<U> coeffs(up_degree()+1,static_cast<U>(0));
        for (size_type i = std::max(deg,0);i<= deg;++i){
            coeffs[i] = operator[](i);
        }
        return mtls::polynomial<U>(coeffs);
    }

    extended_polynomial<U> integrate() const
    {
       extended_polynomial<U> i_data(up_degree()+1,down_degree()+1);
       // Choose integration constant such that P(0) = 0.
       for (size_t i = deg; i <= up_degree(); ++i)
       {
            if(i != -1)
                i_data[i+1] = this->operator[](i)/static_cast<U>(i+1);
            else
                i_data[i] = 0;
       }
       return i_data;
    }

    std::vector<U> const& data() const
    {
        return mtls::polynomial<U>::m_data;
    }
    std::vector<U> & data()
    {
        return mtls::polynomial<U>::m_data;
    }
private:
    template <typename V>
    extended_polynomial<U> &addition(const extended_polynomial<V> & X){

        const int rdd = std::min(down_degree(),X.down_degree());
        const int rud = std::max(up_degree(),X.up_degree());
        extended_polynomial<U> result(rud,rdd);

        int mdd = down_degree();
        if(mdd > rdd){
            for(size_type i=rdd;i<mdd;++i)
                result[i] = X[i];
        }
        else{
            mdd = X.down_degree();
            for(size_type i=rdd;i<mdd;++i)
                result[i] = this->operator[](i);
        }
        int mud = up_degree();
        if(mud < rud){
            for(size_type i = mud+1;i<=rud;++i)
                result[i] = X[i];
        }
        else{
            mud = X.up_degree();
            for(size_type i=mud+1;i<=rud;++i)
                result[i] = this->operator[](i);
        }
        for(size_type i = mdd;i<=mud;++i){
            result[i] = X[i]+operator[](i);
        }
        return (*this = result);
    }

    template <typename V>
    extended_polynomial<U> &multiplication(const V & a){
        for(auto & x:mtls::polynomial<U>::m_data)
            x *= a;
        return *this;
    }

    template <typename V>
    extended_polynomial<U> &division(const V & a){
        for(auto & x:mtls::polynomial<U>::m_data)
            x /= a;
        return *this;
    }
};

template <typename U>
/*returns antipolinomial, the result: F(y) = \int_0^{\infty}{e^(-yx)P(x)dx}*/
auto ExponentIntegral(const mtls::polynomial<U> P){
    extended_polynomial<U> result(P.size(),-P.degree()-1);

    size_t fact = 1;
    for(size_t i=0;i<P.size();++i){
        result[-i-1] = U(fact)*P[i];
        fact *= i;
    }
    return result;
}
#endif
