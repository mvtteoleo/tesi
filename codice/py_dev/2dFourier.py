import numpy as np
import matplotlib.pyplot as plt

N = int(2**3)

x0 = 0
y0 = 0
xN = 2*np.pi
yN = 2*np.pi

x = np.linspace(x0, xN, N, endpoint=False)
y = np.linspace(y0, yN, N, endpoint=False)

X, Y = np.meshgrid(x, y)

"""
lap(u) = f that is solved in Fourier space, then is applied 
the pseudo IBM fix
"""
u = lambda x, y: np.sin(x)*np.sin(y)
f_x = lambda x, y: -u(x, y)
f_y = lambda x, y: -u(x, y)


