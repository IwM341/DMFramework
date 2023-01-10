#ifndef GRID_H
#define GRID_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>

#include "grid_tools.hpp"



namespace Function{

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
class VectorGrid;

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
        h = (N>1 ? (b-a)/(N-1) : 0);
    }

    template <typename U = T>
    UniformGrid(const VectorGrid<U> & VG):a(VG._a()),b(VG._b()),N(VG.size()){
        h = (N>1 ? (b-a)/(N-1) : 0);
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
            return grid.at(i);
        }
    };
    typedef  iterator const_iterator;

    iterator begin() const {return iterator(*this,0);}
    iterator end() const {return iterator(*this,N);}
    iterator last() const {return iterator(*this,N-1);}

    const_iterator cbegin() const {return iterator(*this,0);}
    const_iterator cend() const {return iterator(*this,N);}
    const_iterator clast() const {return iterator(*this,N-1);}

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
        if(N <= 1)
            return 0;
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
    iterator last() const{return (Grid.begin() + (Grid.size()-1));}

    typedef typename std::vector<T>::const_iterator const_iterator;
    auto cbegin() const{return Grid.cbegin();}
    auto cend() const{return Grid.cend();}
    const_iterator clast() const{return (Grid.cbegin() + (Grid.size()-1));}
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

template <typename T>
UniformGrid<T> ApplyRange(const UniformGrid<T> & grid,size_t i0,size_t i1){
    return UniformGrid<T>(grid[i0],grid[i1],i1-i0);
}
template <typename T>
VectorGrid<T> ApplyRange(const VectorGrid<T> & grid,size_t i0,size_t i1){
    return VectorGrid<T>(std::vector<T>(grid.begin() + i0,grid.begin() + i1));
}

};
#endif
