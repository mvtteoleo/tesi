#ifndef HH_FUNCS__HH
#define HH_FUNCS__HH

#include <cmath>
#include <math.h>
#include <vector>

// Forcing term definition
/*
* This equation takes as input x, y, z as x (vector) and t
* and returns the value of the exact solution
*/
template<typename T>
T u_exact (const std::vector<T> &x,const T &t)
{
   return std::sin(x[0])*std::sin(x[1])*std::sin(t);
} 
template<typename T>
T f(const std::vector<T> &x,const T &t)
{
    return (std::cos(t) - 2)*std::sin(x[0])*std::sin(x[1]);
}

#endif // !HH_FUNCS__HH
