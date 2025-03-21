import numpy as np
import matplotlib.pyplot as plt

# Number of points
N = int(1e3)
eps = 0.7

# Define geometrical limits
x0 = -3
xN = 3
h = (xN - x0) / (N - 1)
x = np.linspace(x0, xN, N, endpoint=False) 

# !!! NOTE !!!, this offset may be different from 0, 
# but in that case the IBM is not working anymore as 
# it enforcces the function to be 0 in that region !!
offset=0

# Define functions
u = lambda x: x**3 - eps**2
f_0 = lambda x: 6*x**1 # f = u'' !!!
# u = lambda x: x**2 - eps**2
# f_0 = lambda x: 2 + 0*x
chi = lambda x: 0.5 * (1 - np.sign(abs(x - offset) - eps))
f = lambda x: f_0(x) * (1 - chi(x))
u_ex = lambda x: u(x) * (1 - chi(x))

# Compute exact solution
exact = u_ex(x)

# Laplacian operator with periodic BCs
A_lap = (np.diag(-2*np.ones(N)) + np.diag(np.ones(N-1), -1) + 
         np.diag(np.ones(N-1), 1)) / h**2

# Apply periodic BCs
A_lap[0, -1] = 1/h**2
A_lap[-1, 0] = 1/h**2

# IBM Modification (Immersed Boundary term)
eta = h**2  # Set eta proportional to grid size for stability
chi_vec = chi(x) / eta
B = np.diag(chi_vec)  # IBM term only modifies chi > 0 regions

# Solve (A + B) U = b
A = A_lap + B
b = f(x)

# Impose value on left side
A[0, :] = np.zeros(N)
A[0, 0] = 1
b[0] = u_ex(x0)
# Impose value on right side
A[-1, :] = np.zeros(N)
A[-1, -1] = 1
b[-1] = u_ex(xN)

# Solve the system
Uibm = np.linalg.solve(A, b)

# Plot results
plt.plot(x, Uibm, label='IBM sol')
plt.plot(x, exact, label='Exact sol')
plt.axvline(-eps, color='grey', linestyle="-")
plt.axvline(+eps, color='grey', linestyle="-")
plt.grid()
plt.legend()
plt.show()

