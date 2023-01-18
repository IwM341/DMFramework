#ifndef GRID_OBJECT_H
#define GRID_OBJECT_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>

#include "grid.hpp"
#include "templates.hpp"


namespace Function {

template <typename V,typename GridType,typename DerivedType>
struct GridObject{
    typedef DerivedType Derived;
    GridType Grid;
    std::vector<V> values;

    GridObject(){}
    GridObject(const GridType &Grid,const std::vector<V> &values):Grid(Grid),values(values){}
    GridObject(GridType &&Grid, std::vector<V> &&values):Grid(std::move(Grid)),
        values(std::move(values)){}

    GridObject(const GridType &Grid,size_t N):Grid(Grid),values(N){}
    GridObject(size_t N,const GridType &Grid):Grid(Grid),values(N){}
    GridObject(size_t N,GridType &&Grid):Grid(std::move(Grid)),values(N){}


#define make_operator_gridobject_function(op)   template <typename V1,typename GridType1,typename DerivedType1>\
                                            inline DerivedType operator op(const GridObject<V1,GridType1,DerivedType1> &Y){\
                                    GridObject G(Grid,values.size());\
                                    if(!(Y.Grid == Grid) || Y.values.size()!= values.size()){\
                                        throw std::invalid_argument("different GridObjects in operator");\
                                    }\
                                    else{\
                                    for(size_t i=0;i<values.size();++i)\
                                        G.values[i] = values[i] op Y.values[i];\
                                    return G;\
                                    }\
                                }

    make_operator_gridobject_function(+)
    make_operator_gridobject_function(*)
    make_operator_gridobject_function(/)
    make_operator_gridobject_function(-)


#define make_inc_operator_gridobject_function(op)     template <typename V1,typename GridType1,typename DerivedType1>\
                                        inline DerivedType & operator op(const GridObject<V1,GridType1,DerivedType1> &Y){\
                                    if(!(Y.Grid == Grid) || Y.values.size()!= values.size()){\
                                        throw std::invalid_argument("different GridObjects in operator");\
                                    }\
                                    else{\
                                    for(size_t i=0;i<values.size();++i)\
                                        values[i]  op Y.values[i];\
                                    return *static_cast<DerivedType*>(this);\
                                    }\
                                }

    make_inc_operator_gridobject_function(*=)
    make_inc_operator_gridobject_function(-=)
    make_inc_operator_gridobject_function(+=)
    make_inc_operator_gridobject_function(/=)


    #define make_increment_operator_function(op) template <typename MultType>\
                        inline  DerivedType & operator op(const MultType &val){\
                            for(size_t i=0;i<values.size();++i)\
                                values[i] op val;\
                            return *static_cast<DerivedType*>(this);\
                        }

    make_increment_operator_function(+=)
    make_increment_operator_function(-=)
    make_increment_operator_function(*=)
    make_increment_operator_function(/=)

    #define make_operator_function(op)     template <typename MultType>\
                                inline DerivedType  operator op(const MultType &val) const{\
                                    GridObject G(Grid,values.size());\
                                    for(size_t i=0;i<values.size();++i)\
                                        G.values[i] = values[i] op val;\
                                    return G;\
                                }

    make_operator_function(+)
    make_operator_function(-)
    make_operator_function(*)
    make_operator_function(/)


    inline DerivedType  operator -()const{
        GridObject G(Grid,values.size());
        for(size_t i=0;i<values.size();++i)
            G.values[i] = -values[i];
        return G;
    }

    template <typename OperandType>
    friend inline DerivedType operator +(const OperandType &val,const DerivedType& F){
        return F + val;
    }

    template <typename OperandType>
    friend inline DerivedType operator -(const OperandType &val,const DerivedType& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val - F.values[i];
        return G;
    }

    template <typename OperandType>
    friend inline DerivedType operator *(const OperandType &val,const DerivedType& F){
        return F * val;
    }
    template <typename OperandType>
    friend inline DerivedType operator /(const OperandType &val,const DerivedType& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val / F.values[i];
        return G;
    }

};

template <typename new_type,typename...>
struct GridExtractorType;

template <typename...>
struct GridExtractor;

template <typename V,typename...Other>
struct GridExtractor<GridObject<V,Other...>>{
    template <typename new_value_type>
    using type_template = GridObject<new_value_type,Other...>;
};






template <typename new_type,typename V,typename...Other>
struct GridExtractorType<new_type,GridObject<V,Other...>>{
    typedef GridObject<new_type,Other...> type;
};

template <typename GridType,typename InterpolatorType,typename...>
struct BindGridType;

};

#endif
