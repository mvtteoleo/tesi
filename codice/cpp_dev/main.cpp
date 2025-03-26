#include "mesh.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <vector>
int main (int argc, char *argv[]) {
    // Domain definition
    unsigned int Ndim = 2;
    unsigned int Nx= 5e2;
    unsigned int Ny= 5e2;

    // Iterate to see convergence
    for(double ii = 2; ii<6; ++ii){
        int N = 10*std::pow(2, ii); 
        Nx = N;
        Ny = N;
    // Define x0 and xN, points that define the domain as  
    //  ^ y
    //  | _____. xN
    //  |      |
    //  |      |
    //  |._____| -> x
    // x0
    //
    std::vector<double> x0(Ndim, 0);
    std::vector<double> xN(Ndim, 2*M_PI); //  10);
      

    double t = 0.0;

    // Define the 2D mesh & 2d tensor to store the solution
    Tensor2D<double> mesh(Nx, Ny, 0);

    double h = (xN[1] - x0[1])/Nx;

    Tensor2D<double> final = mesh;
    
    // Domain initialization 
    // Set u_0 on \Omega
    final.applyIC(t, h);
    Tensor2D<double> RK3 = final;

    // Solver
    double dt = 1e-5;//std::pow(10, -ii);
    double Tfinal = 200*dt;
    while(t<Tfinal){
        final.applyBC_ext_dom(t,h);
        final.ExplEuler(dt, h, t);
        RK3.applyBC_ext_dom(t,h);
        RK3.StandardRK3(dt, h, t);
        t+=dt;
    };
        
    ;
    std::cout << "\n Error with N = " << Nx << "| L2 Euler: "<< final.computeL2Error(h, t) <<  " | L2 RK3 : "<< RK3.computeL2Error(h,t) << "\n"; 

        
    }//end for
    
    


    // Export data to python for visualization
    return 0;
}
