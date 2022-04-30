#ifndef PHASE_H
#define PHASE_H

#include <cmath>
#include <string>
class vec3{
public:
	double x,y,z;
	
	inline vec3():x(0),y(0),z(0){}
	
	inline vec3(double x,double y,double z):
		x(x),y(y),z(z) {}
	
	inline static vec3 Polar(double r,double theta,double phi){
		return vec3(r*cos(phi)*sin(theta),r*sin(phi)*sin(theta),r*cos(theta));
	}
	
	inline static vec3 PolarCos(double r,double cosTheta,double phi){
		double sinTheta = sqrt(1-cosTheta*cosTheta);
		return vec3(r*cos(phi)*sinTheta,r*sin(phi)*sinTheta,r*cosTheta);
	}
	
	inline vec3 operator +(const vec3 &second) const{
		return vec3(x+second.x,y+second.y,z+second.z);
	}
	inline vec3 operator -(const vec3 &second) const{
		return vec3(x-second.x,y-second.y,z-second.z);
	}
	
	inline vec3 operator -() const{
		return vec3(-x,-y,-z);
	}
	
	inline vec3 operator +=(const vec3 &second){
		return vec3(x+=second.x,y+=second.y,z+=second.z);
	}
	inline vec3 operator -=(const vec3 & second){
		return vec3(x-=second.x,y-=second.y,z-=second.z);
	}
	
	inline vec3 operator *(double a) const{
		return vec3(x*a,y*a,z*a);
	}
	inline vec3 operator /(double a) const{
		return vec3(x/a,y/a,z/a);
	}
	inline vec3 operator *=(double a){
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	inline vec3 operator /=(double a){
		x/=a;
		y/=a;
		z/=a;
		return *this;
	}
	
	inline double operator*(const vec3 &second) const{
		return x*second.x+y*second.y+z*second.z;
	}
	
	inline double quad() const{
		return x*x+y*y+z*z;
	}
	
	inline double norm() const{
		return sqrt(quad());
	}
	

};


vec3 inline operator *(double a,const vec3 &P){
	return P*a; 
}

class vec4{
	public:
	double t;
	double x,y,z;
	
	inline vec4():t(0),x(0),y(0),z(0){}
	inline vec4(double t,double x,double y,double z):
		t(t),x(x),y(y),z(z) {}
	inline vec4(double t,const vec3 &X): t(t),x(X.x),y(X.y),z(X.z){}
	
	inline vec4(const vec3 &X,double mass):
		t(sqrt(mass*mass + X*X)),x(X.x),y(X.y),z(X.z){}
	
	inline vec3 vecPart() const{
		return vec3(x,y,z);
	}
	
	
	inline vec4 operator +(const vec4 &second) const{
		return vec4(t+second.t,x+second.x,y+second.y,z+second.z);
	}
	inline vec4 operator -(const vec4 &second) const{
		return vec4(t-second.t,x-second.x,y-second.y,z-second.z);
	}
	
	inline vec4 operator -() const{
		return vec4(-t,-x,-y,-z);
	}
	
	inline vec4 operator +=(const vec4 &second){
		return vec4(t+=second.t,x+=second.x,y+=second.y,z+=second.z);
	}
	inline vec4 operator -=(const vec4 &second){
		return vec4(t-=second.t,x-=second.x,y-=second.y,z-=second.z);
	}
	
	inline vec4 operator *(double a) const{
		return vec4(t*a,x*a,y*a,z*a);
	}
	inline vec4 operator /(double a) const{
		return vec4(t/a,x/a,y/a,z/a);
	}
	inline vec4 operator *=(double a){
		t*=a;
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	inline vec4 operator /=(double a){
		t/=a;
		x/=a;
		y/=a;
		z/=a;
		return *this;
	}
	
	inline double operator*(const vec4 &second) const{
		return t*second.t - x*second.x-y*second.y-z*second.z;
	}
	
	inline double quad() const{
		return t*t- x*x - y*y - z*z;
	}
	

	
};

namespace std{
    inline std::string to_string(vec3 V) {
        return std::string("vec3(") + std::to_string(V.x) + "," +
                                    std::to_string(V.y) + "," +
                                    std::to_string(V.z) + ")";
    }
    inline  std::string to_string(vec4 V) {
        return std::string("vec4(") + std::to_string(V.t) + "," +
                                    std::to_string(V.x) + "," +
                                    std::to_string(V.y) + "," +
                                    std::to_string(V.z) + ")";
    }
};

inline vec4 operator *(double a,const vec4 &P){
	return P*a; 
}

inline double E(double m,double p){
	return sqrt(m*m+p*p);
}
inline double P(double m,double E){
	return sqrt(E*E - m*m);
}
#endif
