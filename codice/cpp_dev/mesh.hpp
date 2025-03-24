#ifndef MESH_HPP
#define MESH_HPP


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
                mesh[i][j] = u_exact(x, t);
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
            mesh[i][0] = u_exact(x, t);
        } 
          
        // BC on upper side
        int j=Ny-1;
       for(int i=0; i<Nx; ++i){
            std::vector<T> x{i*h, j*h};
            mesh[i][j] = u_exact(x, t);
        } 
          
        // BC on left side
       for(int j=0; j<Nx; ++j){
            std::vector<T> x{0, j*h};
            mesh[0][j] = u_exact(x, t);
        } 

        // BC on right side
        int i = Nx-1;
       for(int j=0; j<Nx; ++j){
            std::vector<T> x{i*h, j*h};
            mesh[i][j] = u_exact(x, t);
        } 
    }; //end applyBC_ext_dom


    /*
    * Advances the solution of 1 timestep
    */
    Tensor2D timestep(const T& dt, const T& h, const T& t){
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
        return solution;
    }


    /*
    * Obtain the error in L1 norm
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
};


#endif // MESH_HPP
