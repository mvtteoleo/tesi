#include "funcs.hpp"

template<typename T>

// RK3 scheme (scalar to scalar)
T
RK3(double &x,double &t,double dt){
	double a1 = (64.0/120.0);
	double a2 = (50.0/120.0);
	double a3 = (-34.0/120.0);
	double a4 = (90.0/120.0);
	double a5 = (-34.0/120.0);
	// Y2 = u + ... Step
	double Y2 = x + dt*f(t, x)* a1 ;
	// Y3 = Y2 + ... Step
	double tn = (t + 64.0*dt/120.0);
	double Y3 = Y2 + dt*a2*f( tn, Y2) + dt*a3*f( t, x);
	// u_{n+1} = Y3 + ... Step
	double tnn = (t + 80.0*dt/120.0);
	x = Y2 + dt*a5*f( t, x) + dt*a4*f( tnn, Y3);
	return x; 
}
