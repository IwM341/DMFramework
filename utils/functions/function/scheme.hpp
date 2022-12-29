#ifndef SCHEME_H
#define SCHEME_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>


namespace Function{

template <size_t N>
struct Scheme{
    size_t indicies[N];
    double weights[N];


    constexpr static const size_t size = N;

    friend std::ostream& operator <<(std::ostream& os, const Scheme<N>& Interpol){
                size_t prec = os.precision();
                os <<std::setprecision(17);
                os<<"|";
                for (size_t i=0;i<N;i++) {
                    os << Interpol.weights[i]<<"|";
                }
                os << std::endl;
                os<<"|";
                for (size_t i=0;i<N;i++) {
                    os << Interpol.indicies[i]<<"|";
                }
                os << std::endl;
                os << std::setprecision(prec);
                return os;
        }
};


struct SchemeLeft{
template <typename T,template <typename> typename GridType>
    static inline Scheme<1> interpolate(const GridType<T> &Grid,T x){
        if(x<=Grid._a())
            return Scheme<1>({{0},{1}});
        else if(x>=Grid._b())
            return Scheme<1>({{Grid.size()-1},{1}});
        else{
            return Scheme<1>({{Grid.pos(x)},{1}});
        }
    }
};

struct SchemeClosest{
template <typename T,template <typename> typename GridType>
    static inline Scheme<1> interpolate(const GridType<T> &Grid,T x){
        if(x<=Grid._a())
            return Scheme<1>({{0},{1}});
        else if(x>=Grid._b())
            return Scheme<1>({{Grid.size()-1},{1}});
        else{
            size_t i = Grid.pos(x);
            if(x - Grid.at(i) > Grid.at(i+1) - x)
                return Scheme<1>({{i+1},{1}});
            else
                return Scheme<1>({{i},{1}});
        }
    }
};

struct SchemeHisto{
template <typename T,template <typename> typename GridType>
    static inline Scheme<1> interpolate(const GridType<T> &Grid,T x){
        if(x<Grid._a())
            return Scheme<1>({{0},{0}});
        else if(x>Grid._b())
            return Scheme<1>({{0},{0}});
        else{
            size_t i = Grid.pos(x);
            return Scheme<1>({{i},{1}});;
        }
    }
};

struct SchemeRight{
template <typename T,template <typename> typename GridType>
    static inline Scheme<1> interpolate(const GridType<T> &Grid,T x){
        if(x<=Grid._a())
            return Scheme<1>({{0},{1}});
        else if(x>=Grid._b())
            return Scheme<1>({{Grid.size()-1},{1}});
        else{
            return Scheme<1>({{Grid.pos(x) + 1},{1}});
        }
    }
};
struct LinearInterpolator{
template <typename T,template <typename> typename GridType>
    static inline Scheme<2> interpolate(const GridType<T> &Grid,T x){
        if(x<=Grid._a())
            return Scheme<2>({{0,0},{1,0}});
        else if(x>=Grid._b())
            return Scheme<2>({{Grid.size()-1,Grid.size()-1},{1,0}});
        else{
            size_t i = Grid.pos(x);
            T w = (x-Grid.at(i))/(Grid.at(i+1)-Grid.at(i));
            return Scheme<2>({{i,i+1},{1-w,w}});
        }
    }
};
struct LinearExtrapolator{
template <typename T,template <typename> typename GridType>
    static inline Scheme<2> interpolate(const GridType<T> &Grid,T x) {
        size_t i = Grid.pos(x);
        T w = (x-Grid.at(i))/(Grid.at(i+1)-Grid.at(i));
        return Scheme<2>({{i,i+1},{1-w,w}});
    }
};


struct CubicInterpolator{
template <typename T,template <typename> typename GridType>
    static inline Scheme<4> interpolate(const GridType<T> &Grid,T x) {
        if(x<Grid._a())
            return Scheme<4>({{0,0,0,0},{1,0,0,0}});
        else if(x>Grid._b())
            return Scheme<4>({{Grid.size()-1,0,0,0},{1,0,0,0}});
        else{
            size_t i = Grid.pos(x);
            if(i==0) ++i;
            if(i == Grid.size()-2) --i;


            const T x0 = Grid.at(i-1);
            const T x1 = Grid.at(i);
            const T x2 = Grid.at(i+1);
            const T x3 = Grid.at(i+2);

            return Scheme<4>({{i-1,i,i+1,i+2},
                              {(x-x3)*(x-x2)*(x-x1)/((-x1+x0)*(-x3+x0)*(-x2+x0)),
                              -(x-x3)*(x-x2)*(x-x0)/((x1-x3)*(x1-x2)*(-x1+x0)),
                              (x-x3)*(x-x1)*(x-x0)/((x1-x2)*(-x2+x0)*(x2-x3)),
                              -(x-x2)*(x-x1)*(x-x0)/((x2-x3)*(x1-x3)*(-x3+x0))}});
        }
    }
};

};

#endif
