#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <iostream>
#include "funcs.hpp"

template <typename T>
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

    // Print tensor (for debugging)
    void print() const {
        for (const auto& row : mesh) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }


    void applyIC(const T &t, const T &h){
        for (int i=0; i<Nx; ++i) {
            for (int j=0; j<Nx; ++j) {
                std::vector<T> x{i*h, j*h};
                mesh[i][j] = u_exact(x, t);
            }
        }
    };
};


#endif // MESH_HPP
