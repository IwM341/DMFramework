#ifndef HISTOGRAMM_H
#define HISTOGRAMM_H

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
#include "g_function.hpp"
#include "grid_object.hpp"
#include "scheme.hpp"
#include "templates.hpp"

namespace Function{

template <template <typename...> typename ConstructedResultStruct,typename ...Other>
struct HelperType;

template <template <typename...> typename ConstructedResultStruct,typename GridType>
struct HelperType<ConstructedResultStruct,GridType>{
     using ResultType = ConstructedResultStruct<GridType,SchemeHisto>;
};

template <template <typename...> typename ConstructedResultStruct,typename GridType,typename ...Other>
struct HelperType<ConstructedResultStruct,GridType,Other...>{
    using ResultType =
    typename HelperType<BindTemplateLeft<ConstructedResultStruct,GridType,SchemeHisto>::template ResultType,Other...>::ResultType;
};

template <typename V,typename GridType, typename ...Other>
using FuncHistoType = typename HelperType<BindTemplateLeft<GridFunction,V>::template ResultType,GridType,Other...>::ResultType;

template <typename V,typename GridType, typename ...Other>
struct Histogramm;

template <typename V,typename GridType, typename ...Other>
struct Histogramm : public GridObject<Histogramm<V,Other...>,GridType,Histogramm<V,GridType,Other...>>{
    typedef GridObject<Histogramm<V,Other...>,GridType,Histogramm<V,GridType,Other...>> Base;

    Histogramm(){}
    Histogramm(const GridType &Grid,const std::vector<Histogramm<V,Other...>> &values):Base(Grid,values){}
    Histogramm(GridType &&Grid,std::vector<Histogramm<V,Other...>> &&values):Base(std::move(Grid),std::move(values)){}

    Histogramm(const GridType &Grid):Base(Grid.size()-1,Grid){}
    Histogramm(GridType &&Grid):Base(Grid.size()-1,std::move(Grid)){}

    Histogramm(const GridType &Grid,const Other&...OtherGrids):Base(Grid.size()-1,Grid){
        for(size_t i=0;i<Base::values.size();++i){
            Base::values[i] = Histogramm<V,Other...>(OtherGrids...);
        }
    }

    Histogramm(const Base &AbstractGridObject):Base(AbstractGridObject){}
    Histogramm(Base && AbstractGridObject):Base(std::move(AbstractGridObject)){}

    template<typename FuncType>
    inline void map(const FuncType &F){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i].map( FunctorM(0.5*(this->Grid.at(i)+this->Grid.at(i+1)),F)  );
        }
    }



    template <typename T,typename...Other1>
    static Histogramm sameGrid(const Histogramm<T,GridType,Other1...> & Hist){
        Histogramm ret(Hist.Grid);
        for(size_t i=0;i<ret.values.size();++i){
            ret.values[i] = Histogramm<V,Other...>::sameGrid(Hist.values[i]);
        }
        return ret;
    }
    template <typename T,typename Interpol,typename...Other1>
    static Histogramm sameGrid(const GridFunction<T,GridType,Interpol,Other1...> & Func){
        Histogramm ret(ApplyRange(Func.Grid,0,Func.Grid.size()-1));
        for(size_t i=0;i<ret.values.size();++i){
            ret.values[i] = Histogramm<V,Other...>::sameGrid(Func.values[i]);
        }
        return ret;
    }

    auto toFunction() const{
        return GridFunctionCreator2<SchemeHisto>::Create(Base::Grid,Vector(Base::values.size() + 1,[&](size_t i){
                                                            if(i!=Base::values.size())
                                                             return Base::values[i].toFunction()/(Base::Grid[i+1]-Base::Grid[i]);
                                                            else
                                                                return (Base::values[i-1].toFunction()*=0);
                                                         }));
    }



    template <typename T,typename ...Args>
    bool putValue(V value,T x,Args...OtherPos){
        if(x >= Base::Grid._a() && x <Base::Grid._b())
            return Base::values[Base::Grid.pos(x)].putValue(value,OtherPos...);
        return false;
    }

    size_t num_of_element() const{
        size_t summ = 0;
        for(size_t i=0;i<this->values.size();++i){
            summ += this->values[i].num_of_element();
        }
        return summ;
    }

    template <typename IteratorType>
    IteratorType saveIter(IteratorType start,IteratorType end) const{
        for(size_t i=0;i<this->values.size() && start != end;++i){
            start = this->values[i].saveIter(start,end);
        }
        return start;
    }

    std::vector<V> AllValues() const{
        std::vector<V> Vals(this->num_of_element());
        saveIter(Vals.begin(),Vals.end());
        return Vals;
    }

    template <typename IteratorType>
    IteratorType loadIter(IteratorType start,IteratorType end){
        for(size_t i=0;i<this->values.size() && start != end;++i){
            start = this->values[i].loadIter(start,end);
        }
        return start;
    }

    V summ() const {
        V Sum = 0;
        for(size_t i=0;i<this->values.size();++i){
            Sum += this->values[i].summ();
        }
        return Sum;
    }




    std::string gridStr(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->values.size();++i){
            ret += this->values[i].gridStr(prefix +
                                           std::to_string(this->Grid.at(i))+"-"+
                                           std::to_string(this->Grid.at(i+1)) + "\t")+'\n';
        }
        return ret;
    }

    std::string toString(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->values.size();++i){
            ret += this->values[i].toString(prefix +
                                            std::to_string(this->Grid.at(i))+"-"+
                                            std::to_string(this->Grid.at(i+1))  + "\t")+'\n';
        }
        return ret;
    }


    const V &operator[](size_t i) const{
        size_t j=0;
        size_t summ_element = 0;
        size_t summ_element_last = 0;
        for(;summ_element <=i && j<Base::values.size();++j){
            summ_element_last = summ_element;
            summ_element+= Base::values[j].num_of_element();
        }

        return Base::values[j-1].operator [](i-summ_element_last);
    }
    V &operator[](size_t i) {
        size_t j=0;
        size_t summ_element = 0;
        size_t summ_element_last = 0;
        for(;summ_element <=i && j<Base::values.size();++j){
            summ_element_last = summ_element;
            summ_element+= Base::values[j].num_of_element();
        }

        return Base::values[j-1].operator [](i-summ_element_last);
    }
    const auto pos(size_t i) const{
        size_t j=0;
        size_t summ_element = 0;
        size_t summ_element_last = 0;
        for(;summ_element <=i && j<Base::values.size();++j){
            summ_element_last = summ_element;
            summ_element+= Base::values[j].num_of_element();
        }
        return std::tuple_cat(_T(std::make_pair(Base::Grid.at(j-1),Base::Grid.at(j))),Base::values[j-1].pos(i-summ_element_last));
    }


};


template <typename V,typename GridType>
struct Histogramm<V,GridType>: GridObject<V,GridType,Histogramm<V,GridType>>{

    typedef GridObject<V,GridType,Histogramm<V,GridType>> Base;

    Histogramm(){}
    Histogramm(size_t){}
    Histogramm(const GridType &Grid,const std::vector<V> &values):Base(Grid,values){}
    Histogramm(GridType &&Grid,std::vector<V>  &&values):Base(std::move(Grid),std::move(values)){}

    Histogramm(const GridType &Grid):Base(Grid.size()-1,Grid){
        for(size_t i=0;i<Base::values.size();++i){
            Base::values[i] = V(0);
        }
    }
    Histogramm(GridType &&Grid):Base(Grid.size()-1,std::move(Grid)){
        for(size_t i=0;i<Base::values.size();++i){
            Base::values[i] = V(0);
        }

    }

    Histogramm(const Base &AbstractGridObject):Base(AbstractGridObject){}
    Histogramm(Base && AbstractGridObject):Base(std::move(AbstractGridObject)){}

    template <typename T>
    static Histogramm sameGrid(const Histogramm<T,GridType> & Hist){
        return Histogramm(Hist.Grid);
    }

    template <typename T,typename...Other1>
    static Histogramm sameGrid(const GridFunction<T,Other1...> & Func){
        Histogramm ret(ApplyRange(Func.Grid,0,Func.Grid.size()-1));
        return ret;
    }

    template <typename T>
    bool putValue(V value,T x){
        if(x >= Base::Grid._a() && x < Base::Grid._b()){
            Base::values[Base::Grid.pos(x)] += value;
            return true;
        }
        return false;
    }

    template<typename FuncType_T_V>
    inline void map(const FuncType_T_V &f){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = f(0.5*(this->Grid.at(i)+this->Grid.at(i+1)));
        }
    }

    size_t num_of_element() const{
        return this->values.size();
    }

    template <typename IteratorType>
    IteratorType saveIter(IteratorType start,IteratorType end) const{
        for(size_t i=0;i<this->values.size() && start != end;++i){
            *start = this->values[i];
            ++start;
        }
        return start;
    }
    std::vector<V> AllValues() const{
        std::vector<V> Vals(this->num_of_element());
        saveIter(Vals.begin(),Vals.end());
        return Vals;
    }

    template <typename IteratorType>
    IteratorType loadIter(IteratorType start,IteratorType end){
        for(size_t i=0;i<this->values.size() && start != end;++i){
            this->values[i] = *start;
            ++start;
        }
        return start;
    }

    V summ() const {
        V Sum = 0;
        for(size_t i=0;i<this->values.size();++i){
            Sum += this->values[i];
        }
        return Sum;
    }

    auto toFunction() const{
        return GridFunction<V,GridType,SchemeHisto>(Base::Grid,Vector(Base::values.size() + 1,[&](size_t i){
                                                                if( i<Base::values.size())
                                                                    return Base::values[i]/(Base::Grid.at(i+1)-Base::Grid.at(i));
                                                                else
                                                        return V(0);
                                                            }));
    }
    std::string gridStr(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->values.size();++i){
            ret += prefix + std::to_string(this->Grid.at(i))+"-"+std::to_string(this->Grid.at(i+1)) + '\n';
        }
        return ret;
    }

    std::string toString(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->values.size();++i){
            ret += prefix + std::to_string(this->Grid.at(i))+"-"+
                    std::to_string(this->Grid.at(i+1))  +'\t' +
                    std::to_string(this->values[i]) + '\n';
        }
        return ret;
    }
    const auto & operator [](size_t i) const{
        return Base::values[i];
    }
    auto & operator [](size_t i){
        return Base::values[i];
    }
    const auto pos(size_t i) const{
        return std::make_pair(Base::grid[i],Base::grid[i+1]);
    }
};

template <typename V,typename ...Other>
struct GridExtractor<Histogramm<V,Other...>>{
    template <typename new_value_type>
    using type_template = Histogramm<new_value_type,Other...>;
};

template <typename new_type,typename V,typename...Other>
struct GridExtractorType<new_type,Histogramm<V,Other...>>{
    typedef Histogramm<new_type,Other...> type;
};

};

#endif
