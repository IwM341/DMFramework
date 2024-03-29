#ifndef DEBUGDEF_H
#define DEBUGDEF_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <type_traits>
#include <stdio.h>
#include <list>
#include <set>
#include <regex>
#include "functions/vectoroperators.hpp"
#include <mutex>
template <typename U,typename V>
std::ostream & operator << (std::ostream & os, const std::pair<U,V> & P){
    os << "pair(" << P.first << ", " << P.second << ")";
    return os;
}

namespace debugdefs{




    template <typename T>
    std::string to_debug_string(const T &x){
        std::ostringstream os;
        os <<  x;
        return os.str();
    }

    template<typename T, typename U>
    std::string to_debug_string(const std::pair<T,U> &x){
        return to_debug_string(x.first)+"\t"+to_debug_string(x.second);
    }

    template <typename Tuple,int remain>
    struct __tuple_elements_string_list{
        constexpr static int i = std::tuple_size<Tuple>::value - remain;
        static std::string get_list(Tuple const & _Tp,const std::string &delim){
            return to_debug_string(std::get<i>(_Tp)) + delim +
                    __tuple_elements_string_list<Tuple,remain-1>::get_list(_Tp,delim);
        }
    };

    template <class Tuple>
    struct __tuple_elements_string_list<Tuple,0>{
        static std::string get_list(Tuple const & _Tp,const std::string &delim){
            return "";
        }
    };

    template<typename...Args>
    std::string to_debug_string(const std::tuple<Args...> &Tp){
        return "(" +
                __tuple_elements_string_list<
                    std::tuple<Args...>,
                    std::tuple_size<std::tuple<Args...>>::value
                        >::get_list(Tp,", ") + ")";
    }

};
#define SVAR(x) (std::string(#x) + std::string(" = ") + debugdefs::to_debug_string(x))
#define PVAR(x) std::cout << SVAR(x) <<std::endl

#define PDEL() (std::cout << "-----------------------------------------------" <<std::endl)
#define SCOMPARE(x,y) (SVAR(x) + " vs " + SVAR(y))
#define COMPARE(x,y) std::cout << SCOMPARE(x,y) <<std::endl;

template <typename T>
std::string TypeString(){
    std::string fname = __PRETTY_FUNCTION__ ;
    return fname.substr(35,fname.size()-84);
}
#define TypeToString(type) TypeString<type>()

template <typename EnumClass>
struct EnNameStr{
	template <EnumClass em>
	static std::string _str(){
		std::string enum_name = __PRETTY_FUNCTION__ ;
		return  enum_name.substr(69,enum_name.size()-135);
	}
};
#define ValueToString(V) EnNameStr<decltype(V)>::_str<V>()

void print(){
    std::cout <<std::endl;
}
template <typename T>
void print(T data){
    std::cout << data <<std::endl;
}

template <typename ...Args, typename T>
void print(T data,Args...args){
    std::cout << data;
    print(args...);
}

template <typename Delimtype,typename T>
void printd(Delimtype delim,T data){
    std::cout << data <<std::endl;
}

template <typename Delimtype,typename ...Args, typename T>
void printd(Delimtype delim,T data,Args...args){
    std::cout << data << delim;
    printd(delim,args...);
}



class tostream: public std::ofstream{
	std::string _name;
	public:
	const std::string &name() const{return _name;}
	tostream(){
		int Name = rand()*rand() + rand();
		_name = std::to_string(Name) + ".tmp";
		this->open(_name);
	}
	
	template <typename T>
	static tostream save_object(const T & Object){
		tostream stream;
		stream << Object <<std::endl;
		stream.close();
		//std::cout <<"closed" <<std::endl;
		return stream;
	}
	tostream (tostream && stream){
		_name = stream.name();
		stream._name = "";
	}
	
	tostream(const std::string &S){
		_name = S;
		this->open(_name);
	}
	~tostream(){
		if(this->is_open())
			this->close();
		if(_name.size() != 0)
			system((std::string("del ") + _name + "\n").c_str());
	}
};


std::string make_path(const std::string & str){
	return std::string("\"") + std::regex_replace(str, std::regex("\\\\"), "\\\\") + "\"";
}

#ifndef GNUPLOT_PATH
#define GNUPLOT_PATH "gnuplot"
#endif
struct Gnuplot{
	struct PlotData{
		std::string data; // filename or function string
		std::string params;
		PlotData(){}
		PlotData(const std::string &data,const std::string &params):data(data),params(params){}
	};
	std::list<PlotData> plots;
	std::list<tostream> tmp_files;
    std::string show_cmd;
	FILE *gp;
	public:
    Gnuplot(const std::string & path = GNUPLOT_PATH){
        gp = popen(path.c_str(),"w");
		if (!gp) { perror("popen gnuplot"); exit(EXIT_FAILURE); };
        show_cmd = "plot";
	}
	
	template <typename T>
	void plotd(const T & Object,const std::string  & params = ""){
		tmp_files.push_back(tostream::save_object(Object));
		plots.push_back(PlotData(make_path(tmp_files.back().name()),params));
	}
	void plot(const std::string  & body,const std::string  & params = ""){
		plots.push_back(PlotData(body,params));
	}
	void plotf(const std::string & fname, const std::string  & params = ""){
		plots.push_back(PlotData(make_path(fname),params));
	}
	
	void command(const std::string & cmd){
		std::cout << "gnuplot>" << cmd <<std::endl;
		fprintf(gp, "%s\n",cmd.c_str());
		fflush(gp);
	}
	
	
    void show(){
		if(plots.empty() or !gp){
			return;
		}
		auto it = plots.begin();
		
        std::string cmd = show_cmd + " ";

		
		
		
		cmd += (*it).data + " " + (*it).params;
		
		for(++it;it!= plots.end();++it){
			cmd += ", ";
			cmd += (*it).data + " " + (*it).params;
		}
		cmd += "\n";
		std::cout << "gnuplot>" << cmd <<std::endl;
		fprintf(gp, "%s",cmd.c_str());
		
		fflush(gp);
	}
	~Gnuplot(){
		pclose(gp);
	}
};
class pipe_switcher{
    std::map<std::string,std::function<void(const std::string &)>> Funcs;
    std::set<std::string> switch_strings;
    decltype (Funcs.begin()) tmp_func;
    void set_tmp(){
        tmp_func = Funcs.begin();
    }
public:
    pipe_switcher(){
        switch_strings.insert("switch");
        switch_strings.insert("change");
    }

    template <typename T>
    void add_pipe(T * Object, void (T::*Func)(const std::string &) , std::string name){
        if(name == ""){
           name = std::to_string(rand()%10000);
        }
        while(Funcs.find(name) != Funcs.end()){
            name[rand() % (name.size())] = std::to_string(rand()%10)[0];
        }
        Funcs[name] = std::bind(Func,Object,std::placeholders::_1);
        set_tmp();
    }
    template <typename Functype>
    void add_pipe( Functype &&F, std::string name){
        if(name == ""){
           name = std::to_string(rand()%10000);
        }
        while(Funcs.find(name) != Funcs.end()){
            name[rand() % (name.size())] = std::to_string(rand()%10)[0];
        }
        Funcs[name] = F;
        set_tmp();
    }
    std::string pipe_names() const {
        std::stringstream out;
        if(Funcs.empty())
            return "";
        auto it = Funcs.begin();
        out << it->first;
        for(++it;it != Funcs.end();++it)
            out << ", " << it->first;
        return out.str();
    }
    std::string tmp_pipe()const{
        if(!Funcs.empty())
            return tmp_func->first;
        else
            return "no pipes";
    }

    template <typename...Args>
    void set_switch_strings(Args const &...args){
        switch_strings.clear();
        add_switch_strings(args...);
    }
    template <typename...Args>
    inline void add_switch_strings(const std::string &S,Args const &...args){
        switch_strings.insert(S);
        add_switch_strings(args...);
    }
    inline void add_switch_strings(const std::string &S){
        switch_strings.insert(S);
    }

    template <typename...Args>
    inline void remove_switch_strings(const std::string &S,Args const &...args){
        switch_strings.erase(switch_strings.find(S));
        remove_switch_strings(args...);
    }
    inline void remove_switch_strings(const std::string &S){
        switch_strings.erase(switch_strings.find(S));
    }

    bool process_string(const std::string &cmd){
        std::stringstream S(cmd);
        std::string str;
        S >> str;
        if(switch_strings.find(str) == switch_strings.end()){
            return false;
        }
        else{
            S >> str;
            auto it = Funcs.find(str);
            if(it == Funcs.end()){
                std::cout<< "can't find " << str<<std::endl;
            }
            else{
                tmp_func = it;
            }
            return true;
        }
    }
    void command(const std::string & cmd){
        (*tmp_func).second(cmd);
    }
    void exec(){
        std::string S;
        while(1){
            std::cout << tmp_pipe() << ">";
            std::getline(std::cin,S);
            if(S == "quit" or S == "exit"){
                return;
            }
            if(!process_string(S)){
                command(S);
            }
        }
    }
};

struct pipe_menu{

    bool block_exit;

    struct action{
        std::string name;
        std::function<bool(const std::string &)> activation;
        std::function<void(const std::string &)> process;
        template <typename Functype>
        action(const std::string & activation_substring,const Functype& process,const std::string name = ""):
            activation([activation_substring](const std::string & S){
                return S.substr(0,activation_substring.size()) == activation_substring;
            }), process(process),name(name) {}

        template <typename T>
        action(const std::string & activation_substring,T * Object,void (T::*proc)(const std::string &),const std::string name = ""):
            activation([activation_substring](const std::string & S){
                return S.substr(0,activation_substring.size()) == activation_substring;
            }), process(std::bind(proc,Object,std::placeholders::_1)),name(name) {}

        template <typename Functype1,typename Functype2,
                  typename = std::enable_if_t<std::is_invocable_v<Functype1,const std::string &>>>
        action(const Functype1 & activation,const Functype2& process,const std::string name = ""):
            activation(activation), process(process),name(name) {}
        template <typename Functype1,typename T,
                  typename = std::enable_if_t<std::is_invocable_v<Functype1,const std::string &>>>
        action(const Functype1 & activation,T * Object,void (T::*proc)(const std::string &),const std::string name = ""):
            activation(activation), process(std::bind(proc,Object,std::placeholders::_1)),name(name) {}

        template <typename T1,typename Functype2>
        action(T1 * Object1,bool (T1::*proc1)(const std::string &),const Functype2& process,const std::string name = ""):
            activation(std::bind(proc1,Object1,std::placeholders::_1)), process(process),name(name) {}
        template <typename T1,typename T2>
        action(T1 * Object1,bool (T1::*proc1)(const std::string &),
               T2 * Object2,void (T2::*proc2)(const std::string &),const std::string name = ""):
            activation(std::bind(proc1,Object1,std::placeholders::_1)),
            process(std::bind(proc2,Object2,std::placeholders::_1)),name(name) {}
    };
    std::list<action> action_list;
    std::function<void(void)> dialog_action;
    pipe_menu():block_exit(1),dialog_action([](){}){
        add_action(action("exit",this,&pipe_menu::exit,"exit"));
        add_action(action("quit",this,&pipe_menu::exit,"quit"));
    }

    template <typename...Args>
    void add_action(action act){
        action_list.push_back(act);
    }

    template <typename FuncType>
    void add_dialog_action(const FuncType &F){
        dialog_action= F;
    }
    template <typename T>
    void add_dialog_action(T * Object,void (T::*Func)(void)){
        dialog_action = std::bind(Func,Object);
    }

    void remove_action(const std::string & name){
        action_list.remove_if([&name](const action & act){return (act.name == name);});
    }

    void exec(){
        std::string S;
        while(block_exit){
            dialog_action();
            std::getline(std::cin,S);
            for(auto it = action_list.begin();it != action_list.end();++it){
                if(it->activation(S)){
                    it->process(S);
                    break;
                }
            }
        }
    }
    void exit(const std::string & S = ""){
        block_exit = 0;
    }
};

void show_progress(double progress,const std::string&prefix = "",const std::string&postfix = "",int barWidth = 100){
    int pos = progress*barWidth;
    std::cout.flush();
    std::cout << prefix << "[";
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos) std::cout << "=";
        else if (j == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %" << postfix << "\r";
}

template <typename RoutineType>
void show_for_progress(size_t N,const RoutineType & Routine,int barWidth = 70,
                       const std::string&prefix = "",const std::string&postfix = ""){
    int pos = -1;
    show_progress(0.0,prefix,postfix,barWidth);
    for(size_t i =0;i<N;++i){
        Routine(i);
        double progress = double(i+1)/N;
        int n_pos = barWidth * progress;
        if(n_pos != pos){
            pos = n_pos;
            show_progress(progress,prefix,postfix,barWidth);
        }
    }
    std::cout << std::endl;
}

template <typename RoutineType,typename FuncType_M>
void show_for_progress(size_t N,const FuncType_M &M_Func,const RoutineType & Routine,
                       int barWidth = 70,const std::string&prefix = "",const std::string&postfix = ""){
    int pos = 0;
    show_progress(0,prefix,postfix,barWidth);
    size_t task_size = 0;
    for(size_t i=0;i<N;++i){
        task_size += M_Func(i);
    }
    size_t item = 0;
    for(size_t i =0;i<N;++i){
        size_t M = M_Func(i);
        for(size_t j =0;j<M;++j){
            Routine(i,j);
            ++item;
            double progress = (item)/((double)task_size);
            int n_pos = barWidth * progress;
            if(n_pos != pos){
                pos = n_pos;
                show_progress(progress,prefix,postfix,barWidth);
            }

        }
    }
    std::cout << std::endl;
}

class show_prog{
    std::string prefix;
    std::string postfix;
    std::mutex m;
    int barWidth;
    int _pos;
public:
    show_prog(int barWidth = 100,const std::string&prefix = "",
              const std::string&postfix = ""):
        barWidth(barWidth),prefix(prefix),postfix(postfix){
        _pos = 0;
        _show(0,0);
    }
private:
    void _show(int pos,float progress){
        if(pos !=0)
            std::cout.flush();
        std::cout << prefix << "[";
        for (int j = 0; j < barWidth; ++j) {
            if (j < pos) std::cout << "=";
            else if (j == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0 + 0.5) << " %" << postfix << "\r";//std::endl;
    }
public:
    void show(float progress){
        m.lock();
        int pos = progress*barWidth+0.1;
        if(pos > _pos){
            _pos = pos;
            _show(pos,progress);
        }
        m.unlock();
    }
    void end(){
        std::cout << std::endl;
    }
    ~show_prog(){
        end();
    }
};

class progress_bar{
    std::string prefix;
    std::string postfix;
    std::mutex m;
    size_t  cur_prog,full_prog;
    int barWidth;
    int _pos;
public:
    progress_bar(size_t full_prog,size_t barWidth = 100, std::string prefix = "",
                std::string postfix = ""):
        prefix(std::move(prefix)),postfix(std::move(postfix)),
        cur_prog(0),full_prog(full_prog),barWidth(barWidth),_pos(0){
        show();
    }

    void show(){
        std::cout << prefix << "[";
        for (int j = 0; j < barWidth; ++j) {
            if (j < _pos) std::cout << "=";
            else if (j == _pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(cur_prog * 100.0/full_prog + 0.5) << " %" << std::endl;
    }
    void operator ++(){
        m.lock();
        cur_prog++;
        int new_pos = cur_prog*barWidth/full_prog;
        if(new_pos > _pos){
            while(_pos != new_pos){
                _pos++;
            }
            show();
        }
        m.unlock();
    }
    void operator ++(int){
        ++(*this);
    }
};

void wait(){
    std::cin.get();
}

#define debug_return(expr) PVAR((expr)); return expr;

#endif
