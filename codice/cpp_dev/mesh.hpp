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


    /*
    * Advances the solution of 1 timestep
    */
    void
    ExplEuler(const T& dt, const T& h, const T& t){
        Tensor2D<T> solution(Nx, Ny);  // Create a new tensor
        
        for (size_t i = 1; i < Nx - 1; ++i) {
            for (size_t j = 1; j < Ny - 1; ++j) {
                std::vector<T> x{i * h, j * h};

                // Compute Laplacian using finite differences
                T lapU = ((*this)(i-1, j) - 2 * (*this)(i, j) + (*this)(i+1, j) +
                          (*this)(i, j-1) - 2 * (*this)(i, j) + (*this)(i, j+1)) / (h * h);

                // Update solution
                solution(i, j) = (*this)(i, j) + dt * (f(x, t) - lapU);
            }
        }
        (*this) = solution;
    }

    /*
     *   RK3
     *   !! WORST THAT EXPL EULER CHECK BC IMPOSITION !!
     */
void RK3(const T& dt, const T& h, const T& t) {
    Tensor2D<T> Y2(Nx, Ny), Buffer(Nx, Ny), Y3(Nx, Ny), solution(Nx, Ny);
    
    // Apply BC to this to enforce them in Y2
    (*this).applyBC_ext_dom(t + 64.0/120.0*dt, h);
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<T> x{i * h, j * h};
            T lapU = ((*this)(i-1, j) - 2 * (*this)(i, j) + (*this)(i+1, j) +
                      (*this)(i, j-1) - 2 * (*this)(i, j) + (*this)(i, j+1)) / (h * h);
            T f_un = f(x, t) - lapU;
            Y2(i, j) = (*this)(i, j) + (64.0 / 120.0) * dt * f_un;
        }
    }

    
    // Apply BC to Y2 to enforce them in Y3
    Y2.applyBC_ext_dom(t + 80./120.*dt, h);
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<T> x{i * h, j * h};
            T lapU = ((*this)(i-1, j) - 2 * (*this)(i, j) + (*this)(i+1, j) +
                      (*this)(i, j-1) - 2 * (*this)(i, j) + (*this)(i, j+1)) / (h * h);
            T f_un = f(x, t) - lapU;
            Buffer(i, j) = Y2(i, j) - (34.0 / 120.0) * dt * f_un;
        }
    }
    
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<T> x{i * h, j * h};
            T lapU = ((Y2(i-1, j) - 2 * Y2(i, j) + Y2(i+1, j)) +
                      (Y2(i, j-1) - 2 * Y2(i, j) + Y2(i, j+1))) / (h * h);
            T f_Y2 = f(x, t + (64.0 / 120.0) * dt) - lapU;
            Y3(i, j) = Buffer(i, j) + (50.0 / 120.0) * dt * f_Y2;
        }
    }
    
    for (size_t i = 1; i < Nx - 1; ++i) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            std::vector<T> x{i * h, j * h};
            T lapU = ((Y3(i-1, j) - 2 * Y3(i, j) + Y3(i+1, j)) +
                      (Y3(i, j-1) - 2 * Y3(i, j) + Y3(i, j+1))) / (h * h);
            T f_Y3 = f(x, t + (80.0 / 120.0) * dt) - lapU;
            solution(i, j) = Buffer(i, j) + (90.0 / 120.0) * dt * f_Y3;
        }
    }
    
    (*this) = solution;
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
