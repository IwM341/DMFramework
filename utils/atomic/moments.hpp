#ifndef MOMENTS_H
#define MOMENTS_H

#include "../montecarlo.hpp"
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



#endif
