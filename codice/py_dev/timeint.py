import numpy as np
import matplotlib.pyplot as plt

# du/dt = f
f = lambda t: -np.sin(t)
u_ex = lambda t: np.cos(t)

# Third-order Runge-Kutta solver (corrected)
def rk3_solver(x_old, t, dt):
    k1 = dt * f(t)
    k2 = dt * f(t + dt / 2)
    k3 = dt * f(t + dt)
    return x_old + (k1 + 4 * k2 + k3) / 6

# Standard Third-Order Runge-Kutta Method
def rk3_standard(x_old, t, dt):
    k1 = f(t) 
    k2 = f(t + dt/2) * (x_old + dt/2 * k1)
    k3 = f(t + dt) * (x_old - dt * k1 + 2 * dt * k2)
    return x_old + (dt/6) * (k1 + 4 * k2 + k3)

# Third-order Adams-Bashforth Method
def ab3_solver(u_prev, f_prev, dt):
    return u_prev[-1] + (dt / 12) * (23 * f_prev[-1] - 16 * f_prev[-2] + 5 * f_prev[-3])

def compute_error(N, dt):
    tt = np.linspace(0, Tf, N, endpoint=True)
    
    u_rk = np.zeros(N)
    u_ee = np.zeros(N)
    u_rk3 = np.zeros(N)
    u_rk3_std = np.zeros(N)
    u_ab3 = np.zeros(N)
    
    u_rk[0] = u_ex(0)
    u_ee[0] = u_ex(0)
    u_rk3[0] = u_ex(0)
    u_rk3_std[0] = u_ex(0)
    u_ab3[0] = u_ex(0)
    
    f_vals = [f(tt[0])]
    
    for i, t in enumerate(tt[:-1]):
        u_rk[i+1] = rk3_solver(u_rk[i], t, dt)
        u_ee[i+1] = u_ee[i] + dt * f(t)
        u_rk3[i+1] = rk3_solver(u_rk3[i], t, dt)
        u_rk3_std[i+1] = rk3_standard(u_rk3_std[i], t, dt)
        
        f_vals.append(f(t))
        
        if i >= 2:  # Start using AB3 after the first two steps
            u_ab3[i+1] = ab3_solver(u_ab3[i-2:i+1], f_vals[i-2:i+1], dt)
        else:
            u_ab3[i+1] = rk3_solver(u_ab3[i], t, dt)  # Use RK3 for initialization
    
    error_rk = np.linalg.norm(u_rk - u_ex(tt), 2) / np.sqrt(N)
    error_ee = np.linalg.norm(u_ee - u_ex(tt), 2) / np.sqrt(N)
    error_rk3 = np.linalg.norm(u_rk3 - u_ex(tt), 2) / np.sqrt(N)
    error_rk3_std = np.linalg.norm(u_rk3_std - u_ex(tt), 2) / np.sqrt(N)
    error_ab3 = np.linalg.norm(u_ab3 - u_ex(tt), 2) / np.sqrt(N)
    
    return dt, error_rk, error_ee, error_rk3, error_rk3_std, error_ab3

Tf = 2
dt_values = [1e-1, 5e-2, 2.5e-2, 1.25e-2]
errors_rk = []
errors_ee = []
errors_rk3 = []
errors_rk3_std = []
errors_ab3 = []

for dt in dt_values:
    N = int(np.ceil(Tf/dt) + 1)
    dt, error_rk, error_ee, error_rk3, error_rk3_std, error_ab3 = compute_error(N, dt)
    errors_rk.append(error_rk)
    errors_ee.append(error_ee)
    errors_rk3.append(error_rk3)
    errors_rk3_std.append(error_rk3_std)
    errors_ab3.append(error_ab3)

convergence_orders = lambda errors: np.log(np.array(errors[:-1]) / np.array(errors[1:])) / np.log(2)

print("Convergence Order:")
print(f"Euler: {convergence_orders(errors_ee)}")
print(f"RK3: {convergence_orders(errors_rk3)}")
print(f"Standard RK3: {convergence_orders(errors_rk3_std)}")
print(f"AB3: {convergence_orders(errors_ab3)}")

# plt.loglog(dt_values, errors_ee, 's-', label='Error Euler')
plt.loglog(dt_values, errors_rk, 'o-', label='Error Runge Kutta')
plt.loglog(dt_values, errors_rk3, 'd-', label='Error RK3')
# plt.loglog(dt_values, errors_rk3_std, 'x-', label='Error Standard RK3')
# plt.loglog(dt_values, errors_ab3, 'v-', label='Error AB3')

# Reference lines for first, second, and third order convergence
ref_line_1 = [errors_ee[0] * (dt/dt_values[0]) for dt in dt_values]
ref_line_2 = [errors_rk[0] * (dt/dt_values[0])**2 for dt in dt_values]
ref_line_3 = [errors_rk3[0] * (dt/dt_values[0])**3 for dt in dt_values]

# plt.loglog(dt_values, ref_line_1, '--', label='1st Order Ref')
# plt.loglog(dt_values, ref_line_2, '--', label='2nd Order Ref')
plt.loglog(dt_values, ref_line_3, '--', label='3rd Order Ref')

plt.xlabel("Time step dt")
plt.ylabel("L2 norm error")
plt.legend()
plt.grid(True, which="both", linestyle="--")
plt.show()

