#include "mesh.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>
int main (int argc, char *argv[]) {
    // Domain definition
    unsigned int Ndim = 2;
    unsigned int Nx= 1e2;
    unsigned int Ny= 1e2;

    // Iterate to see convergence
    std::cout << "N    |      EE err |     RK3 err"; 
    

    // Control convergence
    std::vector<double> conv1;
    std::vector<double> conv2;


    for(double ii = 0; ii<5; ++ii){
//      int N = 10*std::pow(3, ii); 
//      Nx = N;
//      Ny = N;
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

    Tensor2D<double> test1 = mesh;
    
    // Domain initialization 
    // Set u_0 on \Omega
    test1.applyIC(t, h);
    Tensor2D<double> test2 = test1;

    // Solver
    double Tfinal = 2*M_PI;// 100*dt;
    double dt = 0.1*Tfinal/std::pow(2, ii);// 2e-1
    while(t<Tfinal){

        test1.ExplEuler(dt, h, t);

        test2.StandardRK3(dt, h, t);
        t+=dt;
    };
        
double L2test1 = test1.computeL2Error(h, t);
double L2test2 = test2.computeL2Error(h, t);
    std::cout << std::setprecision(3) << std::scientific 
            << "\n" << dt << 
            "\t"         << L2test1  <<  
            "\t"         << L2test2  << "\n"; 
    conv1.insert(conv1.begin(), L2test1);
    conv2.insert(conv2.begin(), L2test2);
        
    }//end for
    

    // Print Order of convergence
    
    


    // Export data to python for visualization
    return 0;
}
