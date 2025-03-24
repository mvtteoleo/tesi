#include "mesh.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
int main (int argc, char *argv[]) {
    // Domain definition
    unsigned int Ndim = 2;
    unsigned int Nx = 1e1;
    unsigned int Ny = 1e1;

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
      

    double t = 0;

    // Define the 2D mesh // 2d tensor to store the solurtion
//     std::vector<std::vector<double>> tensor(Nx, std::vector<double>(Ny, 0.0));
    Tensor2D<double> mesh(Nx, Ny, 0);
    
    // Domain initialization 
    // Set u_0 on \Omega
    mesh.applyIC(t, h);
    mesh.print();
    // Solver


    // Export data to python for visualization
    return 0;
}
