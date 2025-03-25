#ifndef HH_FUNCS__HH
#define HH_FUNCS__HH

#include <cmath>
#include <math.h>
#include <vector>

/*
* This equation takes as input x, y, z as x (vector) and t
* and returns the value of the exact solution
*/
template<typename T>
T u_exact (const std::vector<T> &x,const T &t)
{
   //return std::sin(x[0])*std::sin(x[1])*std::sin(t);
return ( x[0]*(x[0] - 10) + x[1]*(x[1]-10) )*t;
} 

// Forcing term definition
template<typename T>
T f(const std::vector<T> &x,const T &t)
{
    //return (std::cos(t) - 2)*std::sin(x[0])*std::sin(x[1]);
    return ( x[0]*(x[0] - 1) + x[1]*(x[1]-1) ) + 4*t;
}

#endif // !HH_FUNCS__HH
