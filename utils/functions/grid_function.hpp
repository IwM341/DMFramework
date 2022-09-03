#ifndef GRID_FUNCTION_H
#define GRID_FUNCTION_H

#include "vectoroperators.hpp"
#include "../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>

namespace Function{
template <class T>
class UniformFunction1{
	double a,b;
	size_t N;
	double h;
	std::vector<T> F;
public:

    friend std::ostream& operator <<(std::ostream& os, const UniformFunction1<T>& F){
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.F.size();i++) {
            os << F.a + F.h*i << '\t' <<  debugdefs::to_debug_string(F.F[i])<< std::endl;
		}
		os << std::setprecision(prec);
		return os;
	}


	UniformFunction1(double _a = 0,double _b = 0,size_t _N = 0):a(_a),b(_b),N(_N),F(N){
		if(a>b)
			std::swap(a,b);
		h = (b-a)/(N-1);
	}
	UniformFunction1(size_t _N = 0):a(0),b(1),N(_N),F(_N){
		h = (b-a)/(N-1);
	}
	UniformFunction1(double a,double b,std::vector<T> Y):a(a),b(b),N(Y.size()),F(Y){
		if(a>b){
			std::swap(a,b);
			std::reverse(F.begin(),F.end());
		}
		h = (b-a)/(N-1);
	}
	
	template<typename FuncType>
	UniformFunction1(FuncType f,double _a = 0,double _b = 0,size_t _N = 0):a(_a),b(_b),N(_N),F(N){
		if(a>b)
			std::swap(a,b);
		h = (b-a)/(N-1);
		for(size_t i=0;i<N;++i){
			F[i] = f(a+h*i);
		}
	}
	std::vector<double> Xgrid() const{
		std::vector<double> XG(N);
		for(size_t i=0;i<N;i++){
			XG[i] = a + h*i;
		}
		return XG;
	}
	const std::vector<double> & Fval() const{
		return F;
	}
		
	template<typename FuncType>
	void map(FuncType f){
		for(size_t i=0;i<N;++i)
			F[i] = f(a+i*h);
	}
	
	size_t size()const {return N;}
	
	T operator ()(double x) const{
		int i = (x-a)/h;
		if(x < a)
			return F[0];
		else if(i >= N-1)
			return F[N-1];
		else{
			 double di = (x-a-i*h)/h;
			 return di*F[i+1]+(1-di)*F[i];
		}
	}
};
/*
template <typename T>
class UniformFunctionCubic1 : public UniformFunction1<T>{
	public:
	UniformFunctionCubic1();
	UniformFunctionCubic1(const UniformFunction1 &F):UniformFunction1(F){}
	T operator ()(double x) const{
		int i = (x-a)/h;
		if(x < a)
			return F[0];
		else if(i >= N-1)
			return F[N-1];
		else{
			if(i == 0)
				i = 1;
			else if(i >= N-2)
				i = N-3
				
			double x0 = (x-a-i*h)/h;
			double xm = xi+1.0;
			double xp = 1.0-xi;
			double xpp = 2.0-xi;
			
			return 0.5*xm*xpp*(F[i]*xp+F[i+1]*x0)-x0*xp/6*(F[i-1]*xpp+F[i+2]*xm);
		}
	}
};
*/
extern inline auto find_less(const std::vector<double> &X,double x){
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

template <typename T>
UniformFunction1<double> LoadUniformFunction1(const char *filename,bool isTitle = true){
	std::ifstream ifs(filename, std::ifstream::in);
	std::vector<double> Y;
	double xmin,xmax,x,y;
	if(isTitle){
		std::string S;
		std::getline(ifs,S);
	}
	bool first_time = true;
	while(ifs>>x>>y){
		if(first_time)
			xmin = x;
		xmax = x;
		Y.push_back(y);
	}
    return UniformFunction1<double>(xmin,xmax,Y);
}

template <class T>
class GridFunction1{
	std::vector<double> X;
	std::vector<T> Y;
	
	friend std::ostream& operator<<(std::ostream& os, const GridFunction1<T>& F){
		
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.size();i++) {
            os << F.X[i] << '\t' <<  debugdefs::to_debug_string(F.Y[i])<< std::endl;
		}
		os << std::setprecision(prec);
		return os;
	}
	
public:

	GridFunction1(size_t N = 0):X(N),Y(N){}
	enum SortState {SORTED,FLIPPED,CHAOTIC};
	
private:
	static void PrepareVectors(std::vector<double> &A,std::vector<T> &B,SortState SS){
		if(SS == SORTED)
			return;
		else if(SS == FLIPPED){
			std::reverse(A.begin(), A.end());
			std::reverse(B.begin(), B.end());
		}
		else{
			std::vector<std::pair<double,T>> AB(A.size());
			for(size_t i=0;i<AB.size();i++){
				AB[i] = std::pair<double,T>(A[i],B[i]);
			}
			std::sort(AB.begin(),AB.end(),
				[](const std::pair<double,T> &a,const std::pair<double,T> &b){
					return a.first<b.first;
				}
			);
			for(int i=0;i<AB.size();i++){
				A[i] = AB[i].first;
				B[i] = AB[i].second;
			}
		}
	}
	static void PrepareVector(std::vector<double> &A,SortState SS){
		if(SS == SORTED)
			return;
		else if(SS == FLIPPED)
			std::reverse(A.begin(), A.end());
		else{
			std::sort(A.begin(),A.end());
		}
	}
public:
	GridFunction1(const std::vector<double> &X):X(X),Y(X.size()){
		std::sort(this->X.begin(),this->X.end());
	}
	
	template <typename FuncType>
	GridFunction1(const std::vector<double> &X,FuncType f,bool sorted = false):
	X(X){
		PrepareVector(this->X);
		Y = vmap(f,this->X);
	}
	GridFunction1(const std::vector<double> &X,const std::vector<T> &Y,SortState SS = SORTED):
	X(X),Y(Y){
		PrepareVectors(this->X,this->Y,SS);
	}
	
	double &getX(size_t i){return X[i];}
	double &getY(size_t i){return Y[i];}
	
	const std::vector<double> &getXgrid() const{
		return X;
	}
	
	const std::vector<T> &getYgrid() const{
		return Y;
	}
	
	const size_t size()const{return X.size();}
	T operator ()(double x) const{
		size_t N = X.size();
		size_t i = find_less(X,x);
		size_t i1 = i+1;
		if(i1 >= N){
			return Y[i];
		}
		double a = (x-X[i])/(X[i1]-X[i]);
		return (1-a)*Y[i]+a*Y[i1];
	}
};

template <typename T>
std::vector<T> parse_string(const std::string &S){
	std::stringstream ss(S);
	std::vector<T> parsed;
	T it;
	while(!ss.eof()){
		ss >> it;
		parsed.push_back(it);
	}
	return parsed;
}

extern inline std::map<std::string, std::vector<double>> CSVTable(const std::string & filename){
	
	std::map<std::string, std::vector<double>> Funcs;
	
	std::ifstream ifs(filename, std::ifstream::in);
	
	if(!ifs.is_open()){
		throw std::runtime_error(std::string("CSVTable error: no such file: ") + filename);
	}
	
	std::string S;
	std::getline(ifs,S);
	std::vector<std::string> cols = parse_string<std::string>(S);
	
	for(auto title : cols){
		Funcs[title] = std::vector<double>();
	}
	
	std::vector<double> nums;
	while(!ifs.eof()){
		std::getline(ifs,S);
		nums = parse_string<double>(S);
		if(nums.size()>=cols.size()){
			for(size_t i=0;i<cols.size();i++){
				Funcs[cols[i]].push_back(nums[i]);
			}
		}
	}
	return Funcs;
};

template <class T = double>
class AbstractGrid{
public:
    typedef T value_type;
     size_t size()const;
     T _a()const;
     T _b()const;
     size_t pos(T x) const;
     T at(size_t i) const;
};

template <class T>
class UniformGrid:public AbstractGrid<T>{
    T a;
    T b;
    size_t N;
    T h;
public:
    UniformGrid(T a =0,T b = 1,size_t N = 2):a(a),b(b),N(N){
        h = (b-a)/(N-1);
    }
    size_t size()const {return N;}
    T _a()const {return a;}
    T _b()const {return b;}
    T _h()const {return h;}

    void set_a(T new_a){
        a = new_a;
        h = (b-a)/(N-1);
    }
    void set_b(T new_b){
        b = new_b;
        h = (b-a)/(N-1);
    }
    void set_ab(T new_a,T new_b){
        a = new_a;
        b = new_b;
        h = (b-a)/(N-1);
    }
    void set_N(size_t new_N){
        N = new_N;
        h = (b-a)/(N-1);
    }

    size_t pos(T x) const{
        int i = static_cast<int>( (x-a)/h );
        if(i<0)
            return 0;
        else if(i>= N-1)
            return N-2;
        else
            return i;
    }
    T at(size_t i)const{
        return a + h*i;
    }
};


template <class T = double>
class VectorGrid: public AbstractGrid<T>{

    std::vector<T> Grid;
    void warning(){
        if(Grid.size() < 2)
        std::cout << "warning: Grid.size < 2\n";
    }
public:
    VectorGrid(T a =0,T b = 1,size_t N = 2):
        Grid(Vector(N,[a,b,N](size_t i){return a + i*(b-a)/(N-1);}))
    {
        warning();
    }
    VectorGrid(std::vector<T> Grid):Grid(Grid){
        warning();
    }

    template <typename U>
    VectorGrid(const UniformGrid<U> &unG):Grid(
                                              Vector(unG.size(),[a = unG._a(),b = unG._b(),N =  unG.size()](size_t i){return a + i*(b-a)/(N-1);})
                                              ){}

    size_t size()const {return Grid.size();}
    T _a()const {return Grid[0];}
    T _b()const {return Grid[Grid.size()-1];}


    size_t pos(T x) const{
        size_t i1 = 0;
        size_t i2 = size()-1;
        while(i1 +1 < i2){
            size_t i = (i1 + i2)/2;
            if(x < Grid[i]){
                i2 = i;
            }
            else {
                i1 = i;
            }
        }
        return i1;
    }
    T at(size_t i)const{
        return Grid[i];
    }
};

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


struct LinearInterpolator{
template <typename T,template <typename> typename GridType>
    static inline Scheme<2> interpolate(const GridType<T> &Grid,T x){
        if(x<Grid._a())
            return Scheme<2>({{0,0},{1,0}});
        else if(x>Grid._b())
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

template<typename V,typename GridType,typename InterpolatorType>
struct __GridFunction1{
    GridType Grid;
    std::vector<V> values;

    __GridFunction1(GridType Grid,std::vector<V> values):Grid(Grid),values(values){}
    __GridFunction1(GridType Grid):Grid(Grid),values(Grid.size()){}

    template<typename FuncType_T_V>
    __GridFunction1(GridType Grid,FuncType_T_V F):Grid(Grid),values(Grid.size()){
        for(size_t i =0;i< values.size();++i){
            values[i] = F(Grid.at(i));
        }
    }

    template<typename T>
    V operator ()(T x)const{
        auto Int = InterpolatorType::interpolate(Grid,x);
        V sum = 0;
        for(size_t j = 0;j< Int.size;++j){
            sum += values[Int.indicies[j]]*Int.weights[j];
        }
        return sum;
    }

};



/*GridObject*/
template <typename V,typename GridType>
struct GridObject{
    GridType Grid;
    std::vector<V> values;

    GridObject(){}
    GridObject(GridType Grid,std::vector<V> values):Grid(Grid),values(values){}
    GridObject(GridType Grid,size_t N):Grid(Grid),values(N){}


    #define make_increment_operator_function(op) template <typename T> \
                        inline  GridObject & operator op(const T &val){\
                            for(size_t i=0;i<values.size();++i)\
                                values[i] op val;\
                            return *this;\
                        }

    make_increment_operator_function(+=)
    make_increment_operator_function(-=)
    make_increment_operator_function(*=)
    make_increment_operator_function(/=)

#define make_operator_function(op)     template <typename T>\
                                inline GridObject  operator op(const T &val){\
                                    GridObject G(Grid,values.size());\
                                    for(size_t i=0;i<values.size();++i)\
                                        G.values[i] = values[i] op val;\
                                    return G;\
                                }

    make_operator_function(+)
    make_operator_function(-)
    make_operator_function(*)
    make_operator_function(/)


    inline GridObject  operator -(){
        GridObject G(Grid,values.size());
        for(size_t i=0;i<values.size();++i)
            G.values[i] = -values[i];
        return G;
    }
    template <typename T>
    friend inline GridObject operator +(const T &val,const GridObject& F){
        return F + val;
    }
    template <typename T>
    friend inline GridObject operator -(const T &val,const GridObject& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val - F.values[i];
        return G;
    }
    template <typename T>
    friend inline GridObject operator *(const T &val,const GridObject& F){
        return F * val;
    }
    template <typename T>
    friend inline GridObject operator /(const T &val,const GridObject& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val / F.values[i];
        return G;
    }

#define make_operator_gridobject_function(op)     template <typename T>\
                                inline GridObject  operator op(const GridObject<T,GridType> &Y){\
                                    GridObject G(Grid,values.size());\
                                    if(Y.Grid != Grid || Y.values.size()!= values.size()){\
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


#define make_inc_operator_gridobject_function(op)     template <typename T>\
                                inline GridObject & operator op(const GridObject<T,GridType> &Y){\
                                    if(Y.Grid != Grid || Y.values.size()!= values.size()){\
                                        throw std::invalid_argument("different GridObjects in operator");\
                                    }\
                                    else{\
                                    for(size_t i=0;i<values.size();++i)\
                                        values[i]  op Y.values[i];\
                                    return *this;\
                                    }\
                                }

    make_inc_operator_gridobject_function(*=)
    make_inc_operator_gridobject_function(-=)
    make_inc_operator_gridobject_function(+=)
    make_inc_operator_gridobject_function(/=)
};


/*GridObject*/


template <typename T,typename FuncType>
struct FunctorM{
    T x;
    FuncType F;

    FunctorM(T x,FuncType F):x(x),F(F){}
    template <typename ...Args>
    auto operator ()(Args...args){
        return F(x,args...);
    }


};


template <size_t N,typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction;

template <size_t N,typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction: public GridObject<GridFunction<N-1,V,Other...>,GridType> {
    typedef GridObject<GridFunction<N-1,V,Other...>,GridType> Base;

    GridFunction(){}
    GridFunction(GridType Grid,std::vector<GridFunction<N-1,V,Other...>> values):
        Base(Grid,values){}
    GridFunction(GridType Grid):
        Base(Grid,Grid.size()){}
    GridFunction(const Base & AbstractGridObject):Base(AbstractGridObject){}
    GridFunction(Base && AbstractGridObject) noexcept:Base(std::move_if_noexcept(AbstractGridObject)){}

    template<typename ...GridTypes>
    GridFunction(GridType Grid,GridTypes...OtherGrids):Base(Grid,Grid.size()){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = GridFunction<N-1,V,Other...>(OtherGrids...);
        }
    }

    template <typename FuncType>
    GridFunction(GridType Grid,FuncType F):Base(Grid,Grid.size()){
        for(size_t i =0;i< this->values.size();++i){
            this->values[i] = F(this->Grid.at(i));
        }
    }

    template<typename T,typename ...OtherArgs>
    V operator ()(T x,OtherArgs...args)const{
        auto Int = InterpolatorType::interpolate(this->Grid,x);
        V sum = 0;
        for(size_t j = 0;j< Int.size;++j){
            sum += this->values[Int.indicies[j]](args...)*Int.weights[j];
        }
        return sum;
    }

    template<typename FuncType>
    inline void map(FuncType F){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i].map(  FunctorM(this->Grid.at(i),F)  );
        }
    }

    size_t num_of_element()const{
        size_t summ = 0;
        for(size_t i=0;i<this->values.size();++i){
            summ += this->values.num_of_element();
        }
        return summ;
    }

    std::string gridStr(const std::string & prefix = "") const{
        std::string ret;
        for(size_t i=0;i<this->Grid.size();++i){
            ret += this->values[i].gridStr    (prefix + this->Grid.at(i) + "\t");
        }
        return ret;
    }

};



template <typename V,typename GridType,typename InterpolatorType>
struct GridFunction<1, V, GridType, InterpolatorType> : public GridObject<V,GridType>{
    typedef GridObject<V,GridType> Base;
    InterpolatorType Interpolator;

    GridFunction(){}
    GridFunction(GridType Grid,std::vector<V> values):Base(Grid,values){}
    GridFunction(GridType Grid):Base(Grid,Grid.size()){}

    GridFunction(const Base & AbstractGridObject):Base(AbstractGridObject){}
    GridFunction(Base && AbstractGridObject):Base(std::move_if_noexcept(AbstractGridObject)){}

    template<typename FuncType_T_V>
    GridFunction(GridType Grid,FuncType_T_V F):Base(Grid,Grid.size()){
        for(size_t i =0;i< this->values.size();++i){
            this->values[i] = F(this->Grid.at(i));
        }
    }


    template<typename FuncType_T_V>
    inline void map(FuncType_T_V f){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = f(this->Grid.at(i));
        }
    }

    template<typename T>
    V operator ()(T x)const{
        auto Int = InterpolatorType::interpolate(this->Grid,x);
        V sum = 0;
        for(size_t j = 0;j< Int.size;++j){
            sum += this->values[Int.indicies[j]]*Int.weights[j];
        }
        return sum;
    }

    size_t num_of_element() const{
        return this->values.size();
    }
    std::string gridStr(const std::string & prefix = ""){
        std::string ret;
        for(size_t i=0;i<this->Grid.size();++i){
            ret += this->gridStr(prefix + this->Grid.at(i) + "\t");
        }
        return ret;
    }

};

template <size_t N,typename V,typename GridType, typename ...Other>
struct Histogramm{
    typedef GridObject<Histogramm<N-1,V,Other...>,GridType> Base;

    GridType Grid;
    std::vector<Histogramm<N-1,V,Other...>> values;
    Histogramm(GridType Grid,std::vector<Histogramm<N-1,V,Other...>> values):Grid(Grid),values(values){}
    Histogramm(Base && AbstractGridObject):Base(std::move_if_noexcept(AbstractGridObject)){}

    template <typename T,typename ...Args>
    void putValue(V value,T x,Args...OtherPos){
        if(x >= Grid._a() && x < x._b())
            values[Grid.pos(x)].putValue(value,OtherPos...);
    }

    size_t num_of_element() const{
        size_t summ = 0;
        for(size_t i=0;i<this->values.size();++i){
            summ += this->values.num_of_element();
        }
        return summ;
    }

};


template <typename V,typename GridType>
struct Histogramm<1,V,GridType>{

    typedef GridObject<V,GridType> Base;

    GridType Grid;
    std::vector<V> values;
    Histogramm(GridType Grid):Grid(Grid),values(Grid.size()-1,0){}
    Histogramm(GridType Grid,std::vector<V> values):Grid(Grid),values(values){}
    Histogramm(Base && AbstractGridObject):Base(std::move_if_noexcept(AbstractGridObject)){}


    template <typename T>
    void putValue(V value,T x){
        if(x >= Grid._a() && x < x._b()){
            values[Grid.pos(x)] += value;
        }
    }

    size_t num_of_element() const{
        return this->values.size();
    }
};



};
#endif
