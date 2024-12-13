import numpy as np
import variables_position as vp
import variables_conditions as vc
import Diff
import time
import matplotlib.pyplot as plt

# Start timing
start = time.perf_counter()

# Constants and initial conditions
j, l, k = 0.0, 0.0, 0.0
dp = 0.1 * j
brho = vc.brho0 * (1 + dp / 100)
emi = l * 20
a_max = np.sqrt(emi / vc.betax) * 0.95
a_bin = 2 * a_max / 15
a_init = k * a_bin - a_max

# Derived quantities
dispersion = -3.3494 * dp**6 + 2.2088 * dp**5 + 1.5699 * dp**4 - \
    0.731 * dp**3 - 0.1573 * dp**2 + 0.0121 * dp + 7.0737
MoQ = vc.mass / vc.z
mom = vc.c * vc.z * brho * 1E-6
energy = (np.sqrt(mom**2 + (vc.mass * vc.amu)**2) - (vc.mass * vc.amu)) / vc.mass
beta = ((vc.c * vc.z * brho) / (vc.amu * vc.mass)) / \
    np.sqrt(((vc.c * vc.z * brho) / (vc.amu * vc.mass))**2 + 1)
gamma = 1 / np.sqrt(1 - beta**2)
v = vc.c * beta * 1E-6

# Initialize variables
vx, vy = beta * np.sin(a_init / 1000), beta * np.cos(a_init / 1000)
x, y, t = vc.x0 / 1000, vc.y0 / 1000, 0.0  # Convert mm to m
plot_x, plot_y, plot_vx, plot_vy, plot_t = [], [], [], [], []
stop_flag = 0

# Print initial conditions
print('******* Initial conditions *******')
print(f'Proton number: {vc.z}, Mass: {vc.mass}, MoQ: {MoQ}')
print(f'gamma: {gamma}, beta: {beta}')
print(f'Energy: {energy * 1E6} MeV/u')
print('******* Runge-Kutta parameters *******')
print(f'Initial x: {x * 1000} mm, Initial y: {y * 1000} mm')

# Function for Runge-Kutta step
def runge_kutta_step(x, y, vx, vy, gamma, step):
    k1_x = step * Diff.diff_x(vx)
    k1_y = step * Diff.diff_y(vy)
    k1_vx = step * Diff.diff_vx(x, y, vy, gamma)
    k1_vy = step * Diff.diff_vy(x, y, vx, gamma)

    k2_x = step * Diff.diff_x(vx + k1_vx / 2)
    k2_y = step * Diff.diff_y(vy + k1_vy / 2)
    k2_vx = step * Diff.diff_vx(x + k1_x / 2, y + k1_y / 2, vy + k1_vy / 2, gamma)
    k2_vy = step * Diff.diff_vy(x + k1_x / 2, y + k1_y / 2, vx + k1_vx / 2, gamma)

    k3_x = step * Diff.diff_x(vx + k2_vx / 2)
    k3_y = step * Diff.diff_y(vy + k2_vy / 2)
    k3_vx = step * Diff.diff_vx(x + k2_x / 2, y + k2_y / 2, vy + k2_vy / 2, gamma)
    k3_vy = step * Diff.diff_vy(x + k2_x / 2, y + k2_y / 2, vx + k2_vx / 2, gamma)

    k4_x = step * Diff.diff_x(vx + k3_vx)
    k4_y = step * Diff.diff_y(vy + k3_vy)
    k4_vx = step * Diff.diff_vx(x + k3_x, y + k3_y, vy + k3_vy, gamma)
    k4_vy = step * Diff.diff_vy(x + k3_x, y + k3_y, vx + k3_vx, gamma)

    vx_next = vx + (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx) / 6
    vy_next = vy + (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy) / 6
    x_next = x + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6
    y_next = y + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6
    return x_next, y_next, vx_next, vy_next

# Main loop
while stop_flag < 1:
    step = vc.step_time * vc.c * 1E-9 if (-0.005 < y < 0.005) and (x > 0) else vc.step_time * 100 * vc.c * 1E-9
    x, y, vx, vy = runge_kutta_step(x, y, vx, vy, gamma, step)

    plot_x.append(x * 1000)
    plot_y.append(y * 1000)
    # plot_vx.append(vx)
    # plot_vy.append(vy)
    plot_t.append(t*1E9/vc.c)

    t += step
    if (x > 0) and (y < 0):
        stop_flag = 0.5
    if ((y > 0) and (stop_flag == 0.5)) or (t*1E9/vc.c>500):
        stop_flag = 1

print('******* Revolution time *******')
print(f'{t*1E9/vc.c:.4f} ns')

print('*******   The number of iterations   *******')
print(f'{len(plot_t)} times')

# Plot results
plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)
plt.plot(plot_x, plot_y, label="x-y Trajectory")
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(plot_t, plot_x, label="t-x")
plt.plot(plot_t, plot_y, label="t-y")
plt.legend()

# End timing
end = time.perf_counter()
print('******* Execution time *******')
print(f'{end - start:.2f} seconds')
if t*1E9/vc.c>500:
    print('!!! something strange !!!')
    print('!!! Revolution time is too long !!!')
    print('!!! At first, check step_time in variables_conditions !!!')

plt.show()


