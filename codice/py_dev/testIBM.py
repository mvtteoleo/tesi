import numpy as np
import matplotlib.pyplot as plt

"""
Trying to solve u'' = f on [-1, 1]
f such that : u_ex = x^2-eps^2 & u_ex <0 u_ex=0 
BC: u(-1), u(1), u'(-1), u'(1) from u_ex
Actually this is not IBM, but 'penalization method', 
still could be good for the Poisson equation to impose
the phisical bc//Impose the pressure to be null inside
the body of the poroous material
"""
# Number of points (x_0 to x_{N-1})
N = int(3e2)
eps = 0.5

# !!! NOTE !!!, this offset may be different from 0, 
# but in that case the IBM is not working anymore as 
# it enforcces the function to be 0 in that region !!
offset=0

# Define geometrical limits
x0 = -1
xN = 1
h = (xN-x0)/(N-1)

# Set functions
u    = lambda x: x**2 - eps**2
# f = u''
f_0  = lambda x: 2 + 0*x  # 0*x needed for the plot at least
# chi needs to be 0 in the fluid region and 1 elsewhere
chi  = lambda x: 0.5 * (1 - np.sign(abs(x - offset) - eps))
u_ex = lambda x: u(x) * (1 - chi(x)) 

x = np.linspace(x0, xN, N, endpoint=True)

# Get exact solution in a vector for ease of use
exact = np.zeros(N)
exact = u_ex(x)

# initialize b_ibm
b = np.zeros(N)
b = f_0(x)

# %% Implement IBM
"""
To implement IBM the system moves from u'' = f to 
u'' + chi*u/eta = f that discretized is:
A_lap u + B u = f_b
where B_ij = chi(x_i)* delta_ij /eta
so A U = b
"""

# Laplacian operator
A_lap = (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1) + 
    np.diag(np.ones(N-1), 1) )/h**2

# Needs to be small, if kept as h**2 the  solution is very 
# oscillating probably due to losing the diagonal dominance
eta = h**4
B = np.zeros([N, N])
chi_vec = chi(x)/eta
B = np.diag( np.diag(chi_vec) )

for i in range(2, N-2):
    A_lap[i, i] += B[i]

# IMPOSE BCS

# Dirichlet BCs
A_lap[0, :]  = np.zeros(N)
A_lap[0, 0] = 1
b[0] = u_ex(x0)

A_lap[-1, :]  = np.zeros(N)
A_lap[-1, -1] = 1
b[-1] = u_ex(xN)

Uibm = np.zeros(N)
Uibm = np.linalg.solve(A_lap, b)

# # Print Matrix datas
# print(f"{A_lap = }")
# print(f"{b = }")
# print(f"{Uibm = }")

# Visualize matrix better
# plt.pcolormesh(A_lap, cmap='plasma', shading='auto')
# Modified version of the matrix to better understand the effect of the IBM
plt.pcolormesh(h**2*A_lap, cmap='magma', shading='auto')
plt.colorbar(label="Value")
plt.title("Matrix Visualization using pcolormesh")
plt.show()


# PLOT the result
plot_result=1
if plot_result==1:
    plt.plot(x, Uibm ,'-', label='IBM sol')
    plt.plot(x, exact, label='exact sol')
    plt.axvline(-eps, color='grey')
    plt.axvline(+eps, color='grey')
    plt.grid()
    plt.legend()
    plt.show()

"""Works, the 1d problem without the IBM returns the correct approximation"""
noIBM=1
if noIBM==0:
    # Assemble linear system A*U=b for the simple non IBM case
    A = np.zeros([N, N])
    b = np.zeros(N)
    U = b
    # initialize b
    for i in range(N):
        b[i] = f_0(i*h+x0)

    A = (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1) + np.diag(np.ones(N-1), 1) )/h**2
    # impose BCs

    # Dirichlet
    A[ 0,  :] = np.zeros(N)
    A[ 0,  0] = 1  
    b[0]      = u_ex(x0)
    A[-1,  :] = np.zeros(N)
    A[-1, -1] = 1
    b[-1]     = u_ex(1)

    # Neumann
    A[ 1,  :] = np.zeros(N)
    A[ 1,  0] = -1/h  
    A[ 1,  1] = +1/h
    b[1]      =  u_p(x0)

    A[-2,  :] = np.zeros(N)
    A[-2, -2] = +1/h
    A[-2, -1] = -1/h
    b[-2]     = u_p(xN) 

    # Check impose U(0) = 0
    center0=1
    if center0==0:
        N_2 = int(np.ceil(N/2))
        A[N_2, :] = np.zeros(N)
        A[N_2, N_2] = 1

    # print(A)


    U = np.linalg.solve(A, b)
    # Plot for check
    plot=1
    if plot==0:
        xx = np.linspace(-1, 1, N)
        plt.plot(xx, chi(xx) , label='chi')
        plt.plot(xx, u_ex(xx), label='u_ex')
        plt.plot(xx, f(xx)   , label='f')
        plt.plot(xx, U       , label='U_h')
        plt.legend()
        plt.show()


    plot_Uh=1
    if plot_Uh==0:
        plt.plot(range(N), exact, label='exact solution')
        plt.plot(range(N), U, label='Numerical solution')
        plt.legend()
        plt.show()
