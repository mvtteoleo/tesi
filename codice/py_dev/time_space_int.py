import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x, y, t):
    return np.sin(x) * np.sin(y) * np.sin(t)

def forcing_function(x, y, t):
    return np.sin(x) * np.sin(y) * np.cos(t)

def explicit_euler(dt, T, Nx, Ny):
    x = np.linspace(0, 2 * np.pi, Nx)
    y = np.linspace(0, 2 * np.pi, Ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    u = exact_solution(X, Y, 0)  # Initial condition
    t = 0
    
    while t < T:
        u += dt * forcing_function(X, Y, t)
        t += dt
    
    return u, X, Y

def rk3(dt, T, Nx, Ny):
    x = np.linspace(0, 2 * np.pi, Nx)
    y = np.linspace(0, 2 * np.pi, Ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    u = exact_solution(X, Y, 0)  # Initial condition
    t = 0
    
    while t < T:
        k1 = dt * forcing_function(X, Y, t)
        k2 = dt * forcing_function(X, Y, t + dt / 2)
        k3 = dt * forcing_function(X, Y, t + dt)
        
        u += (k1 + 4 * k2 + k3) / 6
        t += dt
    
    return u, X, Y

def compute_error(numerical, exact):
    return np.sqrt(np.mean((numerical - exact) ** 2))

def convergence_test():
    T = np.pi * 2  # Final time
    Nx, Ny = 40, 40  # Grid size
    dt_values = [0.2, 0.1, 0.05, 0.025, 0.0125]
    
    euler_errors, rk3_errors = [], []
    
    for dt in dt_values:
        u_euler, X, Y = explicit_euler(dt, T, Nx, Ny)
        u_rk3, _, _ = rk3(dt, T, Nx, Ny)
        
        exact = exact_solution(X, Y, T)
        euler_errors.append(compute_error(u_euler, exact))
        rk3_errors.append(compute_error(u_rk3, exact))
    
    plt.loglog(dt_values, euler_errors, 'o-', label='Explicit Euler')
    plt.loglog(dt_values, rk3_errors, 's-', label='RK3')
    plt.xlabel('Time Step (dt)')
    plt.ylabel('L2 Error')
    plt.legend()
    plt.grid(True)
    plt.title('Convergence of Explicit Euler and RK3')
    plt.show()
    
    for i in range(len(dt_values) - 1):
        rate_euler = np.log(euler_errors[i+1] / euler_errors[i]) / np.log(dt_values[i+1] / dt_values[i])
        rate_rk3 = np.log(rk3_errors[i+1] / rk3_errors[i]) / np.log(dt_values[i+1] / dt_values[i])
        print(f'dt={dt_values[i+1]:.5f} -> dt={dt_values[i]:.5f}: Euler Rate = {rate_euler:.2f}, RK3 Rate = {rate_rk3:.2f}')

convergence_test()
