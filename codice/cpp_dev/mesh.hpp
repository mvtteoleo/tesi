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
                (*this)(i, j) = solution(i, j) + dt * (f(x, t) - lapU);
            }
        }
    }


    /*
     * Using the Heun's method (III order):
     */
    void 
    RK3(const T dt, const T h, const T t) {
        // Creo i tensori temporanei
        Tensor2D<T> k1(Nx, Ny, 0), k2(Nx, Ny, 0), k3(Nx, Ny, 0);


        // Obtain k1
        for(int i=1; i<Nx-1; ++i){
            for(int j=1; j<Nx-1; ++j){
                // Get position
                std::vector<T> x{i*h, j*h};

                k1(i, j) = (*this)(i, j) + 
                + dt*( laplacian((*this), i, j, h) + 
                + f(x, t)  );
            }
        }

        // Obtain k2
        for(int i=1; i<Nx-1; ++i){
            for(int j=1; j<Nx-1; ++j){
                // Get position
                std::vector<T> x{i*h, j*h};

                k2(i, j) = (*this)(i, j) + 
                dt*(laplacian(k1, i, j, h) +
                + f(x, t + dt/3.)  );
            }
        }

        // Obtain k3
        for(int i=1; i<Nx-1; ++i){
            for(int j=1; j<Nx-1; ++j){
                // Get position
                std::vector<T> x{i*h, j*h};

                k3(i, j) = (*this)(i, j) + 
                dt*(laplacian(k2, i, j, h) +
                + f(x, t + 2.*dt/3.)  );
            }
        }

        // Obtain (*this)
        for(int i=1; i<Nx-1; ++i){
            for(int j=1; j<Nx-1; ++j){
                // Get position
                (*this)(i, j) = (k1(i, j) + 3*k3(i, j) )/4;
            }
        }

    } // END RK3



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
