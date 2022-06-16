#ifndef DEBUGDEF_H
#define DEBUGDEF_H

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>
namespace debugdefs{
    std::string to_debug_string(double x){
        std::ostringstream os;
        os <<  x;
        return os.str();
    }
    std::string to_debug_string(float x){
        std::ostringstream os;
        os <<  x;
        return os.str();
    }

    template<typename T, typename U>
    std::string to_debug_string(std::pair<T,U> x){
        return to_debug_string(x.first)+"\t"+to_debug_string(x.second);
    }

    template<typename T>
    std::string to_debug_string(T x){
        return std::to_string(x);
    }
};
#define SVAR(x) (std::string(#x) + std::string(" = ") + debugdefs::to_debug_string(x))
#define PVAR(x) std::cout << SVAR(x) <<std::endl

template <typename T>
void print(T data){
    std::cout << data <<std::endl;
}

template <typename ...Args, typename T>
void print(T data,Args...args){
    std::cout << data << "\t";
    print(args...);
}



#endif
