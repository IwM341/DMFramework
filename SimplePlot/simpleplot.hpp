#ifndef SIMPLPLOT_H
#define SIMPLPLOT_H



#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "../grid_function.hpp"

#if defined( __linux__)
    std::string rm = "rm";
#elif defined(_WIN32) || defined(_WIN64)
    std::string rm = "del";
#endif

std::string PrepareString(const std::string &str){
    return "\""+str+"\"";
}

struct DataParam{
    std::vector<double> X,Y;
    std::string legend;
    DataParam(){}
    DataParam(const std::vector<double> &X,
              const std::vector<double> &Y,
              const std::string &legend = ""):X(X),Y(Y),legend(legend){}
};

void DrawPlot(const std::string python,
              const std::string PlotterEXE,
              const std::vector<DataParam> &Data,
              const std::string &yscale = "",
              int dpi = 100,
              const std::string &title = "title",
              const std::string &xlabel = "xlabel",
              const std::string &ylabel = "ylabel"){
    std::string cmd = python + " " + PlotterEXE;
    cmd += " -xl " + PrepareString(xlabel) +
            " -yl " + PrepareString(ylabel) +
            " -tl " + PrepareString(title) +
            " -res " + std::to_string(dpi) + " -sl ";

    for(auto pdata : Data){
        std::ofstream file(PrepareString(pdata.legend));
        file << Function::GridFunction1(pdata.X,pdata.Y) << std::endl;
        file.close();
        cmd += " " + pdata.legend;
    }
    std::cout <<cmd <<std::endl;
    system(cmd.c_str());
    for(auto pdata : Data){
        system( (rm + " " + PrepareString(pdata.legend) ).c_str());
    }

}

#endif
