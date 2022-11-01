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
#include <type_traits>

template <template <typename...>  typename InputClass,typename...Types>
struct BindTemplateLeft{
    template <typename...Args>
    using ResultType = InputClass<Types...,Args...>;
};

template <template <typename...>  typename InputClass,typename...Types>
struct BindTemplateRight{
    template <typename...Args>
    using ResultType = InputClass<Args...,Types...>;
};

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

extern inline std::map<std::string, std::vector<double>> CSVTable(std::istream & stream){
	
	std::map<std::string, std::vector<double>> Funcs;
	
    /*
	std::ifstream ifs(filename, std::ifstream::in);
	
	if(!ifs.is_open()){
		throw std::runtime_error(std::string("CSVTable error: no such file: ") + filename);
	}
    */
	std::string S;
    std::getline(stream,S);
	std::vector<std::string> cols = parse_string<std::string>(S);
	
	for(auto title : cols){
		Funcs[title] = std::vector<double>();
	}
	
	std::vector<double> nums;
    while(!stream.eof()){
        std::getline(stream,S);
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
     size_t pos(const T x) const;
     T at(const size_t i) const;
     inline T operator [](const size_t i) const{
         return this->at(i);
     }

     bool operator ==(const AbstractGrid& OtherGrid) const;
};

template <class T>
class UniformGrid:public AbstractGrid<T>{
    T a;
    T b;
    size_t N;
    T h;

    friend std::ostream & operator << (std::ostream & os,const UniformGrid & VG){
        os << "UniformGrid[";
        os << VG[0];
        for(size_t i=1;i<VG.size();++i){
            os << ",\t" <<VG[i];
        }
        os << "]";
        return os;
    }

public:
    UniformGrid(T a =0,T b = 1,size_t N = 2):a(a),b(b),N(N){
        h = (b-a)/(N-1);
    }

    template <typename ValueType,typename IteratorGrid>
    static UniformGrid fromIter(IteratorGrid startGrid,IteratorGrid endGrid){
        size_t N = 0;
        if(startGrid == endGrid){
            return UniformGrid();
        }
        T _a = *startGrid;
        T _b = _a;
        N++;
        for(;startGrid!=endGrid;++startGrid){
            T _b1 = *startGrid;
            if(_b != _b1){
                _b = _b1;
                N++;
            }
        }
        return UniformGrid(_a,_b,N);
    }
    inline bool operator ==(const UniformGrid& OtherGrid) const{
        return (a == OtherGrid.a) && (b == OtherGrid.b) && (N == OtherGrid.N);
    }

    inline bool operator !=(const UniformGrid& OtherGrid)const{
        return (a != OtherGrid.a) || (b != OtherGrid.b) || (N != OtherGrid.N);
    }


    struct iterator{
        UniformGrid grid;
        size_t i;
        iterator(const UniformGrid &grid,size_t i = 0):i(i),grid(grid){}

        bool operator ==(const iterator & it2){
            return (it2.grid == grid && it2.i==i);
        }
        bool operator !=(const iterator & it2){
            return (it2.grid != grid || it2.i!=i);
        }

        iterator &operator ++(){
            ++i;
            return *this;
        }
        iterator &operator ++(int){
            iterator ret = *this;
            ++i;
            return ret;
        }
        iterator &operator +=(int di){
            i += di;
            return *this;
        }
        iterator &operator -=(int di){
            i -= di;
            return *this;
        }
        iterator operator +(int di){
            iterator ret(this->grid,i+di);
            return ret;
        }
        iterator operator -(int di){
            iterator ret(this->grid,i-di);
            return ret;
        }
        T operator *() const{
            return grid.at[i];
        }
    };
    typedef  iterator const_iterator;

    iterator begin() const {return iterator(*this,0);}
    iterator end() const {return iterator(*this,N);}

    iterator cbegin() const {return iterator(*this,0);}
    iterator cend() const {return iterator(*this,N);}

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

    size_t pos(const T x) const{
        int i = static_cast<int>( (x-a)/h );
        if(i<0)
            return 0;
        else if(i>= N-1)
            return N-2;
        else
            return i;
    }
    T at(const size_t i)const{
        return a + h*i;
    }
    inline T operator [](const size_t i) const{
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

    friend std::ostream & operator << (std::ostream & os,const VectorGrid & VG){
        os << "VectorGrid[";
        os << VG[0];
        for(size_t i=1;i<VG.size();++i){
            os << ",\t" <<VG[i];
        }
        os << "]";
        return os;
    }

public:
    typedef typename std::vector<T>::iterator iterator;
    auto begin() const{return Grid.begin();}
    auto end() const{return Grid.end();}

    typedef typename std::vector<T>::const_iterator const_iterator;
    auto cbegin() const{return Grid.cbegin();}
    auto cend() const{return Grid.cend();}
    const std::vector<T> & grid(){return Grid;}
    VectorGrid(const T a =0,const T b = 1,size_t N = 2):
        Grid(Vector(N,[a,b,N](size_t i){return a + i*(b-a)/(N-1);}))
    {
        warning();
    }
    VectorGrid(const std::vector<T> &Grid):Grid(Grid){
        warning();
    }
    VectorGrid(std::vector<T> &&Grid):Grid(std::move(Grid)){
        warning();
    }

    bool operator == (const VectorGrid & OtherGrid) const{
        return Grid == OtherGrid.Grid;
    }

    bool operator != (const VectorGrid & OtherGrid) const{
        return Grid != OtherGrid.Grid;
    }

    template <typename U>
    VectorGrid(const UniformGrid<U> &unG):Grid(
                                              Vector(unG.size(),[a = unG._a(),b = unG._b(),N =  unG.size()](size_t i){return a + i*(b-a)/(N-1);})
                                              ){}

    template <typename IteratorGrid>
    inline static VectorGrid<T> fromIter(IteratorGrid startGrid,IteratorGrid endGrid){
        size_t N = unique_values_sorted(startGrid,endGrid);
        std::vector<T> grid(N);
        for(size_t i=0;i<N;++i){
            grid[i] = *startGrid;
            startGrid++;
        }
        return VectorGrid<T>(std::move(grid));
    }

    size_t size()const {return Grid.size();}
    T _a()const {return Grid[0];}
    T _b()const {return Grid[Grid.size()-1];}


    size_t pos( const T x) const{
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
    T at(const size_t i)const{
        return Grid[i];
    }
    const T & operator [](const size_t i)  const{
        return Grid[i];
    }

    T & operator [](const size_t i)  {
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



/*GridObject*/
template <template <typename...> typename BaseType,typename T>
struct is_template_base_of_helper{
    template <typename...Args>
    static std::true_type convertable(const volatile BaseType<Args...> *);

    static std::false_type convertable(const volatile void *);

    using value_type = typename decltype(convertable(static_cast<T*>(nullptr)))::type;
};

template <template <typename...> typename BaseType,typename T>
struct is_template_base_of: is_template_base_of_helper<BaseType,T>::value_type{};

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


    #define make_increment_operator_function(op) inline  DerivedType & operator op(const V &val){\
                            for(size_t i=0;i<values.size();++i)\
                                values[i] op val;\
                            return *static_cast<DerivedType*>(this);\
                        }

    make_increment_operator_function(+=)
    make_increment_operator_function(-=)
    make_increment_operator_function(*=)
    make_increment_operator_function(/=)

    #define make_operator_function(op)     inline DerivedType  operator op(const V &val){\
                                    GridObject G(Grid,values.size());\
                                    for(size_t i=0;i<values.size();++i)\
                                        G.values[i] = values[i] op val;\
                                    return G;\
                                }

    make_operator_function(+)
    make_operator_function(-)
    make_operator_function(*)
    make_operator_function(/)


    inline DerivedType  operator -(){
        GridObject G(Grid,values.size());
        for(size_t i=0;i<values.size();++i)
            G.values[i] = -values[i];
        return G;
    }

    friend inline DerivedType operator +(const V &val,const DerivedType& F){
        return F + val;
    }

    friend inline DerivedType operator -(const V &val,const DerivedType& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val - F.values[i];
        return G;
    }

    friend inline DerivedType operator *(const V &val,const DerivedType& F){
        return F * val;
    }
    friend inline DerivedType operator /(const V &val,const DerivedType& F){
        GridObject G(F.Grid, F.values.size());
        for(size_t i=0;i<F.values.size();++i)
            G.values[i] = val / F.values[i];
        return G;
    }

};





/*GridObject*/


template <typename T,typename FuncType>
struct FunctorM{
    const T &x;
    const FuncType &F;

    FunctorM(const T &x,const FuncType &F):x(x),F(F){}
    template <typename ...Args>
    inline auto operator () (Args...args)const{
        return F(x,args...);
    }
};


template <typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction;

template <typename V,typename GridType,typename InterpolatorType, typename ...Other>
struct GridFunction: public GridObject<GridFunction<V,Other...>,GridType,GridFunction<V,GridType,InterpolatorType,Other...>> {
    typedef GridObject<GridFunction<V,Other...>,GridType,GridFunction<V,GridType,InterpolatorType,Other...>> Base;

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
        for(;summ_element <i && j<Base::values.size();++j){
            summ_element_last = summ_element;
            summ_element+= Base::values[j].num_of_element();
        }

        return Base::values[j-1].operator [](i-summ_element_last);
    }
    V &operator[](size_t i) {
        size_t j=0;
        size_t summ_element = 0;
        size_t summ_element_last = 0;
        for(;summ_element <i && j<Base::values.size();++j){
            summ_element_last = summ_element;
            summ_element+= Base::values[j].num_of_element();
        }

        return Base::values[j-1].operator [](i-summ_element_last);
    }

    template<typename FuncType>
    inline void map(const FuncType &F){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i].map( FunctorM(this->Grid.at(i),F)  );
        }
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

    class value_iterator{
        typename decltype(Base::values)::iterator it;
        typename GridFunction<V,Other...>::value_iterator it_deep;

    public:
        bool operator == (const value_iterator & vit1){
            return (it == vit1.it && it_deep == vit1.it_deep);
        }
        bool operator != (const value_iterator & vit1){
            return (it != vit1.it || it_deep != vit1.it_deep);
        }
        value_iterator(decltype(it) it,decltype(it_deep) it_deep):it(it),it_deep(it_deep){}
        value_iterator(decltype(it) it):it(it),it_deep((*it).begin()){}
        inline V & operator *(){
            return *it_deep;
        }

        value_iterator & operator ++(){
            if((++it_deep) == (*it).end()){
                ++it;
                it_deep = (*it).begin();
            }
        }
    };
    value_iterator begin(){
        return value_iterator(Base::values.begin());
    }
    value_iterator end(){
        return value_iterator((decltype (value_iterator::it))Base::values.back(),Base::values.back()->end());
    }


    class const_value_iterator{
        typename decltype(Base::values)::const_iterator it;
        typename GridFunction<V,Other...>::const_value_iterator it_deep;

    public:
        bool operator == (const value_iterator & vit1) const{
            return (it == vit1.it && it_deep == vit1.it_deep);
        }
        bool operator != (const value_iterator & vit1) const{
            return (it != vit1.it || it_deep != vit1.it_deep);
        }
        const_value_iterator(decltype(it) it,decltype(it_deep) it_deep):it(it),it_deep(it_deep){}
        const_value_iterator(decltype(it) it):it(it),it_deep((*it).begin()){}
        inline const V & operator *() const{
            return *it_deep;
        }

        const_value_iterator & operator ++(){
            if((++it_deep) == (*it).end()){
                ++it;
                it_deep = (*it).begin();
            }
        }
    };
    const_value_iterator cbegin() const{
        return const_value_iterator(Base::values.cbegin());
    }
    const_value_iterator cend() const{
        return const_value_iterator((decltype (value_iterator::it))Base::values.back(),Base::values.back()->end());
    }
};



template <typename V,typename GridType,typename InterpolatorType>
struct GridFunction<V, GridType, InterpolatorType> : public GridObject<V,GridType,GridFunction<V, GridType, InterpolatorType>>{
    typedef GridObject<V,GridType,GridFunction<V, GridType, InterpolatorType>> Base;
    InterpolatorType Interpolator;

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

    template<typename FuncType_T_V>
    inline void map(const FuncType_T_V &f){
        for(size_t i=0;i<this->values.size();++i){
            this->values[i] = f(this->Grid.at(i));
        }
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

    typedef typename decltype(Base::values)::iterator value_iterator;
    typedef typename decltype(Base::values)::const_iterator const_value_iterator;

    value_iterator begin(){
        return Base::values.begin();
    }
    value_iterator end(){
        return Base::values.end();
    }

    const_value_iterator cbegin() const{
        return Base::values.cbegin();
    }
    const_value_iterator cend() const{
        return Base::values.cend();
    }

};


template <typename...>
struct GridExtractor;


template <typename new_type,typename...>
struct GridExtractorType;


template <typename new_type,typename V,typename...Other>
struct GridExtractorType<new_type,GridFunction<V,Other...>>{
    typedef GridFunction<new_type,Other...> type;
};

template <typename V,typename...Other>
struct GridExtractor<GridFunction<V,Other...>>{
    template <typename new_value_type>
    using type_template = GridFunction<new_value_type,Other...>;
};

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

template <typename GridType,typename InterpolatorType,typename V, typename...Other>
struct BindGridType<GridType,InterpolatorType,GridFunction<V,Other...>>{
    using ResultType = GridFunction<V,GridType,InterpolatorType,Other...>;
};

template <typename T>
struct VectorType;

template <typename T>
struct VectorType<std::vector<T>>{
    using ResultType = T;
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





template <typename GridType>
GridType diffGrid1(const GridType & grid);

template <typename T>
UniformGrid<T> diffGrid1(const UniformGrid<T> & grid){
    return UniformGrid<T>(grid._a() + grid._h()/2,grid._b() - grid._h()/2,grid.size()-1);
}


template <typename T>
VectorGrid<T> diffGrid1(const VectorGrid<T> & grid){
    VectorGrid<T> retGrid(grid.size()-1);
    for(size_t i=0;i<retGrid.size();++i){
        retGrid[i] = (grid.at(i) + grid.at(i+1))/2;
    }
    return retGrid;
}





template <typename V,typename GridType, typename ...Other>
struct Histogramm;


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


template <typename T>
UniformGrid<T> ApplyRange(const UniformGrid<T> & grid,size_t i0,size_t i1){
    return UniformGrid<T>(grid[i0],grid[i1],i1-i0);
}
template <typename T>
VectorGrid<T> ApplyRange(const VectorGrid<T> & grid,size_t i0,size_t i1){
    return VectorGrid<T>(std::vector<T>(grid.begin() + i0,grid.begin() + i1));
}

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

    class value_iterator{
        typename decltype(Base::values)::iterator it;
        typename Histogramm<V,Other...>::value_iterator it_deep;

    public:
        bool operator == (const value_iterator & vit1){
            return (it == vit1.it && it_deep == vit1.it_deep);
        }
        bool operator != (const value_iterator & vit1){
            return (it != vit1.it || it_deep != vit1.it_deep);
        }
        value_iterator(decltype(it) it,decltype(it_deep) it_deep):it(it),it_deep(it_deep){}
        value_iterator(decltype(it) it):it(it),it_deep((*it).begin()){}
        inline V & operator *(){
            return *it_deep;
        }

        value_iterator & operator ++(){
            if((++it_deep) == (*it).end()){
                ++it;
                it_deep = (*it).begin();
            }
        }
    };
    value_iterator begin(){
        return value_iterator(Base::values.begin());
    }
    value_iterator end(){
        return value_iterator((decltype (value_iterator::it))Base::values.back(),Base::values.back()->end());
    }


    class const_value_iterator{
        typename decltype(Base::values)::const_iterator it;
        typename Histogramm<V,Other...>::const_value_iterator it_deep;

    public:
        bool operator == (const value_iterator & vit1) const{
            return (it == vit1.it && it_deep == vit1.it_deep);
        }
        bool operator != (const value_iterator & vit1) const{
            return (it != vit1.it || it_deep != vit1.it_deep);
        }
        const_value_iterator(decltype(it) it,decltype(it_deep) it_deep):it(it),it_deep(it_deep){}
        const_value_iterator(decltype(it) it):it(it),it_deep((*it).begin()){}
        inline const V & operator *() const{
            return *it_deep;
        }

        const_value_iterator & operator ++(){
            if((++it_deep) == (*it).end()){
                ++it;
                it_deep = (*it).begin();
            }
        }
    };
    const_value_iterator cbegin() const{
        return const_value_iterator(Base::values.cbegin());
    }
    const_value_iterator cend() const{
        return const_value_iterator((decltype (value_iterator::it))Base::values.back(),Base::values.back()->end());
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
