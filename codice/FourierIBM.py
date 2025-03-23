import numpy as np
import matplotlib.pyplot as plt
"""
Trying to solve -u'' = f on [-1, 1]
f such that : u_ex = x^2-eps^2 & u_ex <0 u_ex=0 
BC: u(-1), u(1), u'(-1), u'(1) from u_ex
"""
# Number of points (x_0 to x_{N-1}) 
# For efficiency porpouse is better if they are power
# of 2 due to ffts algorithms
N = int(2**9)
eps = 0.0

# !!! NOTE !!!, this offset may be different from 0, 
# but in that case the IBM is not working anymore as 
# it enforcces the function to be 0 in that region !!
offset= np.pi

# Define geometrical limits
x0 = 0
xN = 2*np.pi
# ******* WARNING ******* 
# Here we do fft so 0 and 2pi are the same 
# value so to impose this we need to have the domain close on
# one of the two end !!
h = (xN-x0)/(N)
x = np.linspace(x0, xN, N, endpoint=False)

eta  = h**5

# Set functions
u    = lambda x: np.sin(x)
# f = u''
f_0  = lambda x: -u(x)

# chi needs to be 0 in the fluid region and 1 elsewhere
# chi  = lambda x: 0.5 * (1 - np.sign(abs(x - offset) - eps))
chi = lambda x: 0.5 * (1 - np.tanh(10 * ((x - offset) - eps)))
u_ex = lambda x: u(x)  * chi(x) 


# Plot chi per check
# plt.plot(x, chi(x))
# plt.show()

# Get exact solution in a vector for ease of use
exact = np.zeros(N)
exact = u_ex(x)

# initialize b
b = np.zeros(N)
b = f_0(x)

"""
Test with fft in 1D
u'' + chi/eta u = f in fourier space becomes
(-k^2 + chi^(k)/eta )*u^(k) = f^(k)
Now the problem is to understand what is going on with chi^(k)
In theory chi^ should be represented as a diagonal matrix and 
not couple the modes, but let's  check
The steps are: 
- Compute the N fourier modes
- Compute f^ and chi^
- Divide and get u^
- ifft(u^) to get the final answer
"""
chi_vals = chi(x)
plt.plot(chi_vals)
b = b*chi_vals
plt.plot(b)
plt.show()

# Compute Fourier frequencies
k = np.fft.fftfreq(N, d=h) * 2 * np.pi  # Frequency domain
k[0] = 1e-10  # Avoid division by zero

# Compute FFT of forcing function and IBM term
F_k = np.fft.fft(b)
B_k = np.zeros(N) # np.fft.fft(chi_vals) / eta  # No matrix needed!



plt.plot(np.real(-k**2 + B_k), label='real pt')
plt.plot(np.imag(-k**2 + B_k), label='imag pt')
plt.legend()
plt.show()

# Solve in Fourier space
U_k = F_k / (-k**2 * h**2 + B_k)  
U_k[0] = 0

# Transform back to real space // In the case of u_ex being sin -like  
# the values are in the real plane
U_fft = np.real(np.fft.ifft(U_k))
U_fft = U_fft*chi_vals

# exit()
# Plot results
plt.figure(figsize=(8, 5))
plt.plot(x, U_fft, '-', label='IBM solution (FFT)')
plt.plot(x, u_ex(x), '--', label='Exact solution', alpha=0.7)
plt.axvline(-eps, color='grey', linestyle="--", label="Interface")
plt.axvline(+eps, color='grey', linestyle="--")
plt.xlabel("x")
plt.ylabel("u(x)", rotation=0)
plt.tight_layout()
plt.legend()
plt.grid()
plt.title("IBM Solution via Fourier Transform")
plt.show()
