#include "mesh.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <vector>
int main (int argc, char *argv[]) {
    // Domain definition
    unsigned int Ndim = 2;
    unsigned int Nx = 4e2;
    unsigned int Ny = 4e2;

    for(double ii = 0; ii<5; ++ii){
        int N = 10*std::pow( 2, ii );
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
    std::vector<double> xN(Ndim, 6.28);
    double h = (xN[1] - x0[1])/Nx;
      

    double t = 0.0;

    // Define the 2D mesh // 2d tensor to store the solurtion
//     std::vector<std::vector<double>> tensor(Nx, std::vector<double>(Ny, 0.0));
    Tensor2D<double> mesh(Nx, Ny, 0);
    
    // Domain initialization 
    // Set u_0 on \Omega
    mesh.applyIC(t, h);

    // Solver
    double dt = 1e-5;
    double Tfinal = 1e-3;
    Tensor2D<double> final = mesh;
    while(t<Tfinal){
        t+=dt;
        final.applyBC_ext_dom(t,h);
        //final = final.timestep(dt, h, t);
        final.timestep(dt, h, t);
    };
    double errorL1 = final.computeL1Error(h, t);
    double errorL2 = final.computeL2Error(h, t);
    std::cout << "\n Error with n  = " << ii << ": "<< errorL1 << "| L2: "<< errorL2<< "\n"; 

    Tensor2D<double> Err = mesh;

    size_t totalPoints = Nx * Ny;

    for (size_t i = 0; i < Nx; ++i) {
        for (size_t j = 0; j < Ny; ++j) {
            std::vector<double> x{i * h, j * h};
            double exactValue = u_exact(x, t);
            Err(i, j) = final(i, j) - exactValue;
        }
    }
    

    }//end for


    // Export data to python for visualization
    return 0;
}
