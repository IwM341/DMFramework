#ifndef CELESTIAL_H
#define CELESTIAL_H

#include "../functions/grid_function.hpp"
#include "../functions/vectoroperators.hpp"
#include <string>
#include <cmath>

#define VESC_SUN 2.056e-3
#define VESC_JUPITER 2.03e-4
#define VESC_EARTH 3.73e-5

static const std::map<std::string, double> ME = {
			{"H",1.0},{"H1",1.0},{"Li",7.0},{"He4",4.0},{"He3",3.0},{"C",12.0},{"C12",12.0},{"C13",13.0},
			{"N14",14.0},{"N15",15.0},{"O",16.0},{"O16",16.0},{"O17",17.0},{"O18",18.0},
			{"Ne",20.0},{"Na",23.0},{"Mg",24.0},{"Al",27.0},{"Si",28.0},{"P",31.0},
			{"S",32.0},{"Cl35",35.0},{"Cl36",36.0},{"Ar",40.0},{"K",39.0},{"Ca",40.0},
			{"Sc",45.0},{"Ti",48.0},{"V",51.0},{"Cr",52.0},
			{"Mn",55.0},{"Fe",56.0},{"Ni",59.0},{"Co",59.0},
			};
static const std::map<std::string, int> QE = {
		{"H",1},{"H1",1},{"He4",2},{"He3",2},{"C12",6},{"C13",6},
		{"N14",7},{"N15",7},{"O",8},{"O16",8},{"O17",8},{"O18",8},
		{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
		{"S",16},{"Cl35",18},{"Cl36",18},{"Ar",18},{"K",19},{"Ca",20},
		{"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Ca",40.0},
		{"Mn",25},{"Fe",26},{"Ni",28}
		};

class BodyModel{
	std::map<std::string, std::vector<double>> BM;
	double _VescMin;
	double _VescMax;
	
	BodyModel(const BodyModel& root);
	BodyModel& operator = (const BodyModel&);
	
	
public:
    BodyModel(double VescOnR,const std::string &LoadFileName):
	_VescMin(VescOnR),
	BM(Function::CSVTable(LoadFileName)){
		if(BM.find("Vesc") == BM.end()){
			double phi0 = BM.at("phi")[0];
			double phi1 = BM.at("phi").back();
			_VescMax = _VescMin*sqrt(phi0/phi1);
            BM["Vesc"] = vmap([VescR=_VescMin,phi1](double x)->double
							{
                                return sqrt(x/phi1)*VescR;
							}, BM.at("phi"));
		}
	}
	const std::vector<double> &operator[](const std::string & column) const{
		return BM.at(column);
	}
    double VescMin() const{
		return _VescMin;
	}
	double VescMax() const{
		return _VescMax;
	}
	Function::UniformFunction1<double> UniformRadFunc
						(const std::string & column) const{
							
		return Function::UniformFunction1<double>(0,1,BM.at(column));
	}
};

#endif
