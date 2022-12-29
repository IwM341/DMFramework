#ifndef GRID_TOOLS_H
#define GRID_TOOLS_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>

namespace Function {


template <typename T>
extern inline auto find_less(const std::vector<T> &X,T x){
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
};

#endif
