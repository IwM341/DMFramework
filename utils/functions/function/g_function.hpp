#ifndef G_FUNCTION_H
#define G_FUNCTION_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>

#include "grid_object.hpp"


namespace Function {

template <typename T,typename FuncType>
struct FunctorM{
    const T &x;
    const FuncType &F;

    FunctorM(const T &x,const FuncType &F):x(x),F(F){}
    template <typename ...Args>
    inline auto operator () (Args&&...args)const{
        return F(x,args...);
    }
};


template <typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction;

template <typename new_type,typename V,typename...Other>
struct GridExtractorType<new_type,GridFunction<V,Other...>>{
    typedef GridFunction<new_type,Other...> type;
};

template <typename V,typename...Other>
struct GridExtractor<GridFunction<V,Other...>>{
    template <typename new_value_type>
    using type_template = GridFunction<new_value_type,Other...>;
};


template <typename GridFunctionType>
struct GridFunctionReturnType;

template <typename V,typename ...Other>
struct GridFunctionReturnType<GridFunction<V,Other...>>{
    typedef V type;
};

template <typename LambdaType,typename GridFunctionType>
struct GridLambdaApplyResultType;

template <typename LambdaType,typename V,typename GridType,typename InterpolatorType>
struct GridLambdaApplyResultType<LambdaType,GridFunction<V,GridType,InterpolatorType>>{
    typedef typename std::invoke_result<LambdaType(typename GridType::value_type)>::type type;
};

template <typename LambdaType,typename V,typename GridType,typename InterpolatorType,typename ...Other>
struct GridLambdaApplyResultType<LambdaType,GridFunction<V,GridType,InterpolatorType,Other...>>{
    typedef typename GridLambdaApplyResultType<FunctorM<typename GridType::value_type,LambdaType>,GridFunction<V,Other...>>::type type;
};

template <typename LambdaType,typename GridFunctionType>
struct GridIndexLambdaApplyResultType;

template <typename LambdaType,typename V,typename GridType,typename InterpolatorType>
struct GridIndexLambdaApplyResultType<LambdaType,GridFunction<V,GridType,InterpolatorType>>{
    typedef typename std::invoke_result<LambdaType(size_t)>::type type;
};

template <typename LambdaType,typename V,typename GridType,typename InterpolatorType,typename ...Other>
struct GridIndexLambdaApplyResultType<LambdaType,GridFunction<V,GridType,InterpolatorType,Other...>>{
    typedef typename GridIndexLambdaApplyResultType<FunctorM<size_t,LambdaType>,GridFunction<V,Other...>>::type type;
};



template <typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction: public GridObject<GridFunction<V,Other...>,GridType,GridFunction<V,GridType,InterpolatorType,Other...>> {
    typedef GridObject<GridFunction<V,Other...>,GridType,GridFunction<V,GridType,InterpolatorType,Other...>> Base;
    constexpr static size_t Ndim = GridFunction<V,Other...>::Ndim + 1;
    GridFunction(){}
    GridFunction(const GridType &Grid,const std::vector<GridFunction<V,Other...>> &values):
        Base(Grid,values){}

    GridFunction(GridType &&Grid,std::vector<GridFunction<V,Other...>> &&values):
        Base(std::move(Grid),std::move(values)){}

    GridFunction(const GridType &Grid):
        Base(Grid,Grid.size()){}

    GridFunction(GridType &&Grid):
        Base(Grid.size(),Grid){}

    GridFunction(const Base & AbstractGridObject):Base(AbstractGridObject){}
    GridFunction(Base && AbstractGridObject):Base(std::move(AbstractGridObject)){}


    template<typename ...GridTypes>
    GridFunction(const GridType &Grid,const GridTypes&...OtherGrids):Base(Grid,Grid.size()){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = GridFunction<V,Other...>(OtherGrids...);
        }
    }

    template <typename FuncType>
    GridFunction(const GridType &Grid,const FuncType &F):Base(Grid,Grid.size()){
        for(size_t i =0;i< this->values.size();++i){
            this->values[i] = F(this->Grid.at(i));
        }
    }

    template <typename T>
    static GridFunction sameGrid(const GridFunction<T,GridType,InterpolatorType,Other...> & GFun){
        GridFunction ret(GFun.Grid);
        for(size_t i=0;i<ret.values.size();++i){
            ret.values[i] =GridFunction<V,Other...>::sameGrid(GFun.values[i]);
        }
        return ret;
    }


    //static GridFunction fromString(const std::string)

    template<typename T,typename ...OtherArgs>
    V operator ()(T x,OtherArgs...args)const{
        auto Int = InterpolatorType::interpolate(this->Grid,x);
        V sum = this->values[Int.indicies[0]](args...)*Int.weights[0];
        for(size_t j = 1;j< Int.size;++j){
            sum += this->values[Int.indicies[j]](args...)*Int.weights[j];
        }
        return sum;
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
        return std::tuple_cat(_T(Base::Grid.at(j-1)),Base::values[j-1].pos(i-summ_element_last));
    }

    template<typename FuncType>
    inline void map(const FuncType &F){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i].map( FunctorM(this->Grid.at(i),F)  );
        }
    }

    template<typename FuncType>
    inline void index_map(const FuncType &F){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i].index_map( FunctorM(i,F) );
        }
    }

    template <typename LambdaType>
    auto ApplyLambda(const LambdaType & func){
        using retType = typename GridLambdaApplyResultType<LambdaType,GridFunction>::type;
        auto Ret = GridExtractorType<retType,GridFunction>::type::sameGrid(*this);
        Ret.map(func);
        return Ret;
    }

    template <typename LambdaType>
    auto ApplyIndexLambda(const LambdaType & func){
        using retType = typename GridIndexLambdaApplyResultType<LambdaType,GridFunction>::type;
        auto Ret = GridExtractorType<retType,GridFunction>::type::sameGrid(*this);
        Ret.index_map(func);
        return Ret;
    }
    template <typename LambdaType>
    auto Composition(const LambdaType & Func) const{
        using retType = typename std::result_of<LambdaType(V)>::type;
        typename GridExtractorType<retType,GridFunction>::type Ret(Base::Grid);
        for(size_t i=0;i<Base::values.size();++i){
            Ret.values[i] = Base::values[i].Composition(Func);
        }
        return Ret;
    }

    size_t num_of_element()const{
        size_t summ = 0;
        for(size_t i=0;i<this->values.size();++i){
            summ += this->values[i].num_of_element();
        }
        return summ;
    }

    std::string gridStr(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->Grid.size();++i){
            ret += this->values[i].gridStr(prefix + std::to_string(this->Grid.at(i)) + "\t")+'\n';
        }
        return ret;
    }

    std::string toString(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->Grid.size();++i){

            ret += this->values[i].toString(prefix + std::to_string(this->Grid.at(i)) + "\t")+'\n';
        }

        return ret;
    }

    template <typename IteratorType>
    IteratorType saveIter(IteratorType start,IteratorType end) const{
        for(size_t i=0;i<this->values.size() && start != end;++i){
            start = this->values[i].saveIter(start,end);
        }
        return start;
    }

    template <typename IteratorType>
    IteratorType loadIter(IteratorType start,IteratorType end){
        for(size_t i=0;i<this->values.size() && start != end;++i){
            start = this->values[i].loadIter(start,end);
        }
        return start;
    }

    std::vector<V> AllValues() const{
        std::vector<V> Vals(this->num_of_element());
        saveIter(Vals.begin(),Vals.end());
        return Vals;
    }

    template <typename ...Args>
    static GridFunction fromCSV(size_t start,size_t end,std::vector<double> grid,Args...args){
        GridFunction ret(GridType::fromIter(grid.begin()+start,grid.begin()+end));
        if(start == end){
            return ret;
        }

        size_t js = start;
        size_t jtmp = start;
        size_t i = 0;

        double point_tmp  = grid[jtmp];

        for(;jtmp < end;++jtmp){
            if(point_tmp != grid[jtmp]){
                ret.values[i] = GridFunction<V,Other...>::fromCSV(js,jtmp,args...);
                js = jtmp;
                i++;
            }
        }
        return ret;
    }

    friend V dot(const GridFunction & F,const GridFunction & G){
        V summ = 0;
        for(size_t i=0;i<F.values.size();++i){
            summ += dot(F.values[i],G.values[i]);
        }
        return summ;
    }

};




template <typename V,typename GridType,typename InterpolatorType>
struct GridFunction<V, GridType, InterpolatorType> : public GridObject<V,GridType,GridFunction<V, GridType, InterpolatorType>>{
    typedef GridObject<V,GridType,GridFunction<V, GridType, InterpolatorType>> Base;
    constexpr static size_t Ndim = 1;

    GridFunction(){}
    GridFunction(const GridType &Grid,const std::vector<V> &values):Base(Grid,values){}

    GridFunction(GridType &&Grid,std::vector<V> &&values):Base(std::move(Grid),std::move(values)){}

    GridFunction(const GridType &Grid):Base(Grid,Grid.size()){}
    GridFunction(GridType &&Grid):Base(Grid.size(),std::move(Grid)){}

    GridFunction(const Base & AbstractGridObject):Base(AbstractGridObject){}
    GridFunction(Base && AbstractGridObject):Base(std::move(AbstractGridObject)){}

    template<typename FuncType_T_V>
    GridFunction(const GridType &Grid,const FuncType_T_V &F):
        Base(Grid,Grid.size()){
        for(size_t i =0;i< this->values.size();++i){
            this->values[i] = F(this->Grid.at(i));
        }
    }

    template <typename T>
    static GridFunction sameGrid(const GridFunction<T,GridType,InterpolatorType> & GFun){
        return GridFunction(GFun.Grid);
    }

    const V &operator[](size_t i) const{
        return Base::values[i];
    }
    V &operator[](size_t i) {
        return Base::values[i];
    }
    const auto pos(size_t i) const{
        return _T(Base::Grid.at(i));
    }

    template<typename FuncType_T_V>
    inline void map(const FuncType_T_V &f){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = f(this->Grid.at(i));
        }
    }

    template<typename FuncType_T_V>
    inline void index_map(const FuncType_T_V &f){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = f(i);
        }
    }

    template <typename LambdaType>
    auto ApplyLambda(const LambdaType & func){
        using retType = typename GridLambdaApplyResultType<LambdaType,GridFunction>::type;
        auto Ret = GridExtractorType<retType,GridFunction>::type::sameGrid(*this);
        Ret.map(func);
        return Ret;
    }

    template <typename LambdaType>
    auto ApplyIndexLambda(const LambdaType & func){
        using retType = typename GridIndexLambdaApplyResultType<LambdaType,GridFunction>::type;
        auto Ret = GridExtractorType<retType,GridFunction>::type::sameGrid(*this);
        Ret.index_map(func);
        return Ret;
    }

    template <typename LambdaType>
    auto Composition(const LambdaType & Func) const{
        using retType = typename std::result_of<LambdaType(V)>::type;
        return typename GridExtractorType<retType,GridFunction>::type (Base::Grid,vmap(Func,Base::values));
    }

    V operator ()(typename GridType::value_type x)const{
        auto Int = InterpolatorType::interpolate(this->Grid,x);
        V sum = this->values[Int.indicies[0]]*Int.weights[0];
        for(size_t j = 1;j< Int.size;++j){
            sum += this->values[Int.indicies[j]]*Int.weights[j];
        }
        return sum;
    }

    size_t num_of_element() const{
        return this->values.size();
    }
    std::string gridStr(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->Grid.size();++i){
            ret += prefix + std::to_string(this->Grid.at(i)) + '\n';
        }
        return ret;
    }

    std::string toString(const std::string & prefix = "") const{
        std::stringstream ret;
        for(size_t i=0;i<this->Grid.size();++i){
            ret << prefix << this->Grid.at(i) << '\t' <<  this->values[i] << '\n';
        }
        return ret.str();
    }

    template <typename IteratorType>
    IteratorType saveIter(IteratorType start,IteratorType end) const{
        for(size_t i=0;i<this->values.size() && start != end;++i){
            *start = this->values[i];
            ++start;
        }
        return start;
    }

    template <typename IteratorType>
    IteratorType loadIter(IteratorType start,IteratorType end){
        for(size_t i=0;i<this->values.size() && start != end;++i){
            this->values[i] = *start;
            ++start;
        }
        return start;
    }

    std::vector<V> AllValues() const{
        std::vector<V> Vals(this->num_of_element());
        saveIter(Vals.begin(),Vals.end());
        return Vals;
    }

    static GridFunction fromCSV(size_t start,size_t end,std::vector<double> grid,std::vector<double> values){
        return GridFunction(GridType::fromIter(grid.begin() + start,grid.begin() + end),
                            std::vector<V>(values.begin()+start,values.begin()+end));
    }

    static GridFunction fromCSV(size_t start,size_t end,std::vector<double> grid){
        return GridFunction(GridType::fromIter(grid.begin() + start,grid.begin() + end));
    }

    friend V dot(const GridFunction & F,const GridFunction & G){
        V summ = 0;
        for(size_t i=0;i<F.values.size();++i){
            summ += F.values[i]*G.values[i];
        }
        return summ;
    }
};

template <typename GridType,typename InterpolatorType,typename V, typename...Other>
struct BindGridType<GridType,InterpolatorType,GridFunction<V,Other...>>{
    using ResultType = GridFunction<V,GridType,InterpolatorType,Other...>;
};

template <typename InterpolatorType>
struct GridFunctionCreator2{

    template <typename GridType,typename InternalGridFunctionType>
    static auto Create(const GridType & Grid,const std::vector<InternalGridFunctionType> & values){
        return typename BindGridType<std::remove_reference_t<GridType>,InterpolatorType,InternalGridFunctionType>::ResultType(Grid,values);
    }
    template <typename GridType,typename InternalGridFunctionType>
    static auto Create(GridType && Grid,std::vector<InternalGridFunctionType> && values){
        return typename BindGridType<std::remove_reference_t<GridType>,InterpolatorType,InternalGridFunctionType>::ResultType(std::move(Grid),std::move(values));
    }

    template <typename GridType,typename InitFuncType>
    static auto Create(const GridType & Grid,const InitFuncType& Func){
        return typename BindGridType<std::remove_reference_t<GridType>,InterpolatorType,typename std::result_of_t<InitFuncType(typename std::remove_reference_t<GridType>::value_type)>>::ResultType(Grid,Func);
    }

    template <typename GridType,typename InitFuncType>
    static auto Create(GridType && Grid,const InitFuncType & Func){
        return typename BindGridType<std::remove_reference_t<GridType>,InterpolatorType,typename std::result_of_t<InitFuncType(typename std::remove_reference_t<GridType>::value_type)>>::ResultType(std::move(Grid),Func);
    }
};


template <typename InterpolatorType>
struct GridFunctionCreator1{
    template <typename GridType,typename V>
    static auto Create(const GridType & Grid,const std::vector<V> & values){
        return GridFunction<V,std::remove_reference_t<GridType>,InterpolatorType>(Grid,values);
    }
    template <typename GridType,typename V>
    static auto Create(GridType && Grid,std::vector<V> && values){
        return GridFunction<V,std::remove_reference_t<GridType>,InterpolatorType>(std::move(Grid),std::move(values));
    }

    template <typename GridType,typename InitFuncType>
    static auto Create(const GridType & Grid,const InitFuncType & Func){
        return GridFunction<typename std::result_of<InitFuncType(typename std::remove_reference_t<GridType>::value_type)>::type,std::remove_reference_t<GridType>,InterpolatorType>(Grid,Func);
    }
    template <typename GridType,typename InitFuncType>
    static auto Create( GridType && Grid,const InitFuncType & Func){
        return GridFunction<typename std::result_of<InitFuncType(typename std::remove_reference_t<GridType>::value_type)>::type,std::remove_reference_t<GridType>,InterpolatorType>(std::move(Grid),Func);
    }
};




};

#endif
