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
    for(double ii = 0; ii<6; ++ii){
        int N = 10*std::pow(2, ii );
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
    std::vector<double> xN(Ndim, 8*M_PI); //  10);
    double h = (xN[1] - x0[1])/Nx;
      

    double t = 0.0;

    // Define the 2D mesh & 2d tensor to store the solution
    Tensor2D<double> mesh(Nx, Ny, 0);
    Tensor2D<double> final = mesh;
    
    // Domain initialization 
    // Set u_0 on \Omega
    final.applyIC(t, h);

    // Solver
    double dt = 2e-3;
    double Tfinal = 10*dt;
    while(t<Tfinal){
        final.applyBC_ext_dom(t,h);
        //final.ExplEuler(dt, h, t);
        final.RK3(dt, h, t);
        t+=dt;
    };
        
    double errorL1 = final.computeL1Error(h, t);
    double errorL2 = final.computeL2Error(h, t);
    std::cout << "\n Error with N = " << Nx << ": "<< errorL1 << "| L2: "<< errorL2 << "\n"; 

        
    }//end for
    
    


    // Export data to python for visualization
    return 0;
}
