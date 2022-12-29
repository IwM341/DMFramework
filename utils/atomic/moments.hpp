#ifndef MOMENTS_H
#define MOMENTS_H

#include "../montecarlo/montecarlo.hpp"
#include "../complex/complex_ex.hpp"
#include "../debug/debugdef.hpp"
#include <cmath>

#define Iion 0.2833594103
#define Imgd 0.7165877840
#define Rd 13.6e-9 //GeV

template <typename UniformGenerator01>
MC::MCResult<double> PhiGenerate(UniformGenerator01 G){
	double phi = sqrt(1/sqrt(sqrt(G())) - 1);
	
	double rd = 0;
	if(phi != 0)
		rd = 32.0*exp(-4*atan(phi)/phi)/(3.0*(1-exp(-2*M_PI/phi)));
	
	return MC::MCResult<double>(phi,rd);
}


const std::vector<double> Prob40 = {0.0,.7744047131, .8985894394, 
	.9417436282, .9620050791, .9732002677, .9800581036, .9845710997, 
	.9877023423, .9899652655, .9916545565, .9929494135, .9939640018, 
	.9947739013, .9954307867, .9959709725, .9964205903, .9967988375, 
	.9971200759, .9973952243, .9976327049, .9978390987, .9980196086, 
	.9981783917, .9983188027, .9984435722, .9985549417, .9986547644, 
	.9987445842, .9988256948, .9988991875, .9989659876, .9990268842, 
	.9990825531, .9991335763, .9991804568, .9992236316, .9992634812, 
	.9993003386, .9993344962,1.0};

template <typename UniformGenerator01>
MC::MCResult<int> NGenerate(UniformGenerator01 G){
	double x = G();
	double rd = Imgd;
	
	size_t imn = 2;
	size_t imx = Prob40.size()-1;
	
	while(imn+1<imx){
		size_t i = (imn+imx)/2;
		if(Prob40[i] < x)
			imn = i;
		else
			imx = i;
	}
	
	return MC::MCResult<int>(imn,rd);
}

extern inline double deltaEIon(double phi){
	return Rd*(phi*phi+1.0);
}
extern inline double deltaEMgd(int n){
	return Rd*(1.0-1.0/(n*n));
}

extern inline double phiMax(double deltaEmax){
    return sqrt(deltaEmax/Rd-1.0);
}

extern inline int nMax(double deltaEmax){
    int res = 1.0/sqrt(1.0-deltaEmax/Rd);
    if(res > 100 || deltaEmax>Rd){
        return 100;
    }
    else return res;
}

enum IonFactorDiff{
    dphi,dE_Rd
};
double IonFactor(double s,double phi, IonFactorDiff F = dphi){
    double sum = phi+s;
    double delta = s-phi;
    double diff_fact = (F == dphi ? 2*phi : 1.0);

    double exps = (phi <= 1e-3 ?
                       exp(-4/(1.0+s*s)) :
                       exp(-2*(atan(sum)-atan(delta))/phi)/
                                       (1.0-exp(-2*M_PI/phi)));


    return diff_fact*(128.0*s*s*(1.0+3*s*s+phi*phi))*
                      exps/(3*fast_pow((1+sum*sum)*(1+delta*delta),3));
}

double MigdalFactor(double s, int n){
    return s*s/3*256*fast_pow(n,8)*fast_pow(1.0/(n*n-1.0),4)/(n*(n*n-1.0))*fast_pow((n-1.0)/(n+1.0),2*n);
}
double ElasticFactor(double s){
    return 1.0/fast_pow(1+s*s/4,4);
}


const double _me = 0.52e-3;
const double _mp = 0.938;
const double _alpha = 0.0073;

enum ScatteringType{
    ELASTIC,IONIZATION,MIGDAL
};
enum Target{
    ELECTRON,PROTON
};

template <ScatteringType ST,Target T,typename GenType = void>
struct dF_H{
    MC::MCResult<double> EnergyLoss(double Emax);
    double ScatterFactor(double q);
};

template <Target T,typename GenType>
struct dF_H<MIGDAL,T,GenType>{
    struct atomic_state{
        size_t nMgd;
        double Emgd;

        atomic_state():atomic_state(1){}
        atomic_state(size_t nMgd):nMgd(nMgd),Emgd(deltaEMgd(nMgd)){}
        operator double() const{
            return Emgd;
        }
    };
    GenType &G;
    MC::MCResult<atomic_state> EnergyLoss(double Emax){
        double dnMx = nMax(Emax)-2;
        if(dnMx > 0){
            double nMigdal = 2 + G()*(int)(dnMx+1);
            //deltaE = deltaEMgd(nMigdal);
            return atomic_state(nMigdal,dnMx+1);
        }
        else{
           return 0;
        }
    }
    double ScatterFactor(double q,atomic_state S);
};

template<Target T,typename GenType>
double dF_H<MIGDAL,T,GenType>::ScatterFactor(double q,atomic_state S){

    double s = q/(_alpha*(T == PROTON ? _mp : _me));
    return MigdalFactor(s,S.nMgd);
}

template <Target T>
struct dF_H<ELASTIC,T,void>{
    MC::MCResult<double> EnergyLoss(double Emax){
        return MC::MCResult<double>(0.0,1);
    }
    double ScatterFactor(double q,double EnLoss);
};

template<Target T>
double dF_H<ELASTIC,T,void>::ScatterFactor(double q,double enLoss){
    double s = q/(_alpha*(T == PROTON ? _mp : _me));
    return ElasticFactor(s);
}


template <Target T,typename GenType>
struct dF_H<IONIZATION,T,GenType>{
    GenType &G;
    dF_H(GenType &G):G(G){}
    MC::MCResult<double> EnergyLoss(double Emax){
        if(Emax > Rd){
            //dE = max energy loss by ionization - min energyloss od ionization
            double dE = Emax-Rd;
            return MC::MCResult<double>(Rd + G()*dE,dE/Rd);
        }
        else{
            MC::MCResult<double>(0,0);
        }
    }
    double ScatterFactor(double q,double EnLoss);
};

template<Target T,typename GenType>
double dF_H<IONIZATION,T,GenType>::ScatterFactor(double q,double EnLoss){
    double s = q/(_alpha*(T == PROTON ? _mp : _me));
    if(EnLoss > Rd)
        return IonFactor(s,phiMax(EnLoss),dE_Rd);
    else return 0.0;
}

#endif
