#ifndef MESH_HPP
#define MESH_HPP


#include <cmath>
#include <vector>
#include <iostream>
#include "funcs.hpp"

template <typename T>
/*
* Class to handle the IC, BC and general
* interfaces between solution
*/
class Tensor2D {
private:
    std::vector<std::vector<T>> mesh;
    size_t Nx, Ny;

public:
    // Constructor
    Tensor2D(size_t Nx, size_t Ny, T init_val = T()) 
        : Nx(Nx), Ny(Ny), mesh(Nx, std::vector<T>(Ny, init_val)) {}

    // Access element with (i, j)
    T& operator()(size_t i, size_t j) { return mesh[i][j]; }
    const T& operator()(size_t i, size_t j) const { return mesh[i][j]; }

    // Get dimensions
    size_t get_Nx() const { return Nx; }
    size_t get_Ny() const { return Ny; }

    /* 
     * Print tensor (for debugging)
    */
    void print() const {
        for (const auto& row : mesh) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }


    /*
    *Kinda self-explanatory, takes u_exact 
    *from func.hpp and sets the initial
    * values accordingly
    */
    void applyIC(const T &t, const T &h){
        for (int i=0; i<Nx; ++i) {
            for (int j=0; j<Nx; ++j) {
                std::vector<T> x{i*h, j*h};
                 (*this)(i, j)= u_exact(x, t);
            }
        }
    };

    /*
    * Applyes the boundary conditions on the 
    * sides of the domain  
    */
    void applyBC_ext_dom(const T &t, const T &h){
        // BC on lower side
       for(int i=0; i<Nx; ++i){
            std::vector<T> x{i*h, 0};
            (*this)(i, 0) = u_exact(x, t);
        } 
          
        // BC on upper side
        int j=Ny-1;
       for(int i=0; i<Nx; ++i){
            std::vector<T> x{i*h, j*h};
            (*this)(i, j) = u_exact(x, t);
        } 
          
        // BC on left side
       for(int j=0; j<Nx; ++j){
            std::vector<T> x{0, j*h};
            (*this)(0, j)= u_exact(x, t);
        } 

        // BC on right side
        int i = Nx-1;
       for(int j=0; j<Nx; ++j){
            std::vector<T> x{i*h, j*h};
            (*this)(i, j) = u_exact(x, t);
        } 
    }; //end applyBC_ext_dom

    T
    laplacian(const Tensor2D<T> &u, int i, int j, T h){
        return (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4 * u(i, j) ) / (h * h);
    }

    /*
    * Advances the solution of 1 timestep
    * Using Expicit Euler method
    */
    void
    ExplEuler(const T& dt, const T& h, const T& t){
        Tensor2D<T> solution = (*this);  // Create a new tensor
        solution.applyBC_ext_dom(t, h);
        
        for (size_t i = 1; i < Nx - 1; ++i) {
            for (size_t j = 1; j < Ny - 1; ++j) {
                std::vector<T> x{i * h, j * h};

                // Compute Laplacian using finite differences
                T lapU = laplacian( solution, i, j, h);

                // Update solution
               (*this)(i, j) = solution(i, j) + dt * (f(x, t) - 0*lapU);
            }
        }
    }


    void 
    StandardRK3(const double dt, const double h, const double t) {
    Tensor2D<double> k1(Nx, Ny, 0), k2(Nx, Ny, 0), k3(Nx, Ny, 0), u_temp(Nx, Ny, 0);

    k1.applyBC_ext_dom(t, h);
    // Compute k1 = f(x, y, t^n) - laplacian(u^n)
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<double> x{i * h, j * h};  // Position vector (x, y)
            
            double lapU = ((*this)(i+1, j) + (*this)(i-1, j) +
                          (*this)(i, j+1) + (*this)(i, j-1) -
                          4 * (*this)(i, j)) / (h * h);
            
            k1(i, j) = f(x, t) - 0*lapU;
        }
    }

    // Compute k2 = f(x, y, t^n + dt/2) - laplacian(u^n + dt/2 * k1)
    k2.applyBC_ext_dom(t + 0.5*dt, h);
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<double> x{i * h, j * h};
            
            u_temp(i, j) = (*this)(i, j) + (dt / 2.0) * k1(i, j);
            
            double lapU = 0*(u_temp(i+1, j) + u_temp(i-1, j) +
                           u_temp(i, j+1) + u_temp(i, j-1) -
                           4 * u_temp(i, j)) / (h * h);
            
            k2(i, j) = f(x, t + dt/2) - lapU;
        }
    }

    // Compute k3 = f(x, y, t^n + dt) - laplacian(u^n - dt*k1 + 2*dt*k2)
    k3.applyBC_ext_dom(t + dt, h);
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<double> x{i * h, j * h};
            
            u_temp(i, j) = (*this)(i, j) - dt * k1(i, j) + 2 * dt * k2(i, j);
            
            double lapU = 0*(u_temp(i+1, j) + u_temp(i-1, j) +
                           u_temp(i, j+1) + u_temp(i, j-1) -
                           4 * u_temp(i, j)) / (h * h);
            
            k3(i, j) = f(x, t + dt) - lapU;
        }
    }

    // Update the solution using RK3 formula
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            (*this)(i, j) += (dt / 6.0) * (k1(i, j) + 4.0 * k2(i, j) + k3(i, j));
        }
    }
}



    /*
    * Obtain the error in L1 norm
    * err = sum_i( abs(u_i - u_exact) ) / (N_pts)
    */
    T computeL1Error(T & h, T & t) const {
        T errorSum = 0.0;
        size_t totalPoints = Nx * Ny;

        for (size_t i = 0; i < Nx; ++i) {
            for (size_t j = 0; j < Ny; ++j) {
                std::vector<T> x{i * h, j * h};
                T exactValue = u_exact(x, t);
                T error = std::abs( (*this)(i, j) - exactValue );
                errorSum += error;
            }
        }

        return errorSum / totalPoints; // Normalize by number of points
    }

    /*
    * Obtain the error in L2 norm
    * err = sqrt(  sum( (u_i -u_exact)^2 ) ) / (N_pts)
    */
    T 
    computeL2Error(T & h, T & t) const {
        T errorSum = 0.0;
        size_t totalPoints = Nx * Ny;

        for (size_t i = 0; i < Nx; ++i) {
            for (size_t j = 0; j < Ny; ++j) {
                std::vector<T> x{i * h, j * h};
                T exactValue = u_exact(x, t);
                T error = std::pow( ((*this)(i, j)  - exactValue ), 2);
                errorSum += error;
            }
        }

        return std::sqrt(errorSum )/ totalPoints; // Normalize by number of points
    }
};


#endif // MESH_HPP
