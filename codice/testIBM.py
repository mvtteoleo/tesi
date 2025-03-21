import numpy as np
import matplotlib.pyplot as plt

"""
Trying to solve -u'' = f on [-1, 1]
f such that : u_ex = x^2-eps^2 & u_ex <0 u_ex=0 
BC: u(-1), u(1), u'(-1), u'(1) from u_ex
"""
# Number of points (x_0 to x_{N-1})
N = int(1e2)
eps = 0.5

# Define geometrical limits
a = -1
b = 1
h = (b-a)/(N-1)

# Set functions
u    = lambda x: -x**2 - eps**2
f_0  = lambda x: -2 + 0*x  # 0*x needed for the plot at least
# chi needs to be 0 in the fluid region and 1 elsewhere
chi  = lambda x: 0.5 * (1 - np.sign(abs(x) - eps))
f    = lambda x: f_0(x) * chi(x)
u_ex = lambda x: u(x) * (1 - chi(x)) 
u_p  = lambda x: 2*x

x = np.linspace(a, b, N)

# Plot chi per check
# plt.plot(x, chi(x))
# plt.show()

# Get exact solution in a vector for ease of use
exact = np.zeros(N)
for i in range(N):
    exact[i] = u_ex(i*h + a)


# %% Implement IBM
"""
To implement IBM the system moves from -u'' = f to 
u'' + chi*u/eta = f that discretized is:
A_lap u + B u = f_b
where B_ij = chi(x_i)* delta_ij /eta
so A U = b
"""
# Laplacian operator
A_lap = -(np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1) + 
    np.diag(np.ones(N-1), 1) )/h**2
b = np.zeros(N)

eta = 1/h**4
B = np.zeros([N, N])
chi_vec = chi(x)/eta
B = np.diag( np.diag(chi_vec) )
print(f"{B = }")

# initialize b_ibm
for i in range(N):
    b[i] = f_0(i*h+a)

for i in range(2, N-2):
    A_lap[i, i] -= B[i]

# impose BCs

# Dirichlet periodic
A_lap[0, -1] = -1/h**2
# Neumann
A_lap[-1, -1]= -1/h**2

print(A_lap)
Uibm = np.zeros(N)
Uibm = np.linalg.solve(A_lap, b)
plt.plot(range(N), Uibm , label='IBM sol')
plt.plot(range(N), exact, label='exact sol')
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
        b[i] = f_0(i*h+a)

    A = (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1) + np.diag(np.ones(N-1), 1) )/h**2
    # impose BCs

    # Dirichlet
    A[ 0,  :] = np.zeros(N)
    A[ 0,  0] = 1  
    b[0]      = u_ex(a)
    A[-1,  :] = np.zeros(N)
    A[-1, -1] = 1
    b[-1]     = u_ex(1)

    # Neumann
    A[ 1,  :] = np.zeros(N)
    A[ 1,  0] = -1/h  
    A[ 1,  1] = +1/h
    b[1]      =  u_p(a)

    A[-2,  :] = np.zeros(N)
    A[-2, -2] = +1/h
    A[-2, -1] = -1/h
    b[-2]     = u_p(b) 

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
