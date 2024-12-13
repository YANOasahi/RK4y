import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time

import Bz

# start of an execution timer
start = time.perf_counter()

# *******   physics constant   *******
# light speed [m/s]
c = 299792458
# amu [ev/c^2]
amu = 931494061

# *********************************************
# **************   modify here   **************
# *********************************************
# target nuclide
z, mass = 32.0, 77.9229  # 78Ge
# brho of particles when dp/p=0 [Tm]
brho0 = 4.7447
# beta of the ring in the X-axis
betax = 7.817
# the end time of Runge-Kutta
stop_time = 380  # in ns
# step time of Runge-Kutta
step_time = 0.01  # in ns
# *********************************************

# *******   positions of particles   *******
x0 = 9287.959673
y0 = 0.0
# initial position
r0 = np.array([x0/1000, y0/1000, 0])

# *******   velocities of particles   *******
# for searching an ideal condition
j, l, k = 0.0, 0.0, 0.0
dp = 0.1 * j
brho = brho0 * (1 + dp / 100)
emi = l * 20
a_max = np.sqrt(emi / betax) * 0.95
a_bin = 2 * a_max / 15
a_init = k * a_bin - a_max

# q/m
MoQ = mass/z
# momentum based on brho
mom = c*z*brho*1E-6
# beam energy based on momentum
energy = (np.sqrt(mom**2+(mass*amu)**2)-(mass*amu))/mass
# beta based on brho and q/m
beta = ((c*z*brho)/(amu*mass))/np.sqrt(((c*z*brho) /
                                        (amu*mass))*((c*z*brho)/(amu*mass))+1)
# Lorentz factor
gamma = 1/np.sqrt(1-beta**2)
beta = ((c * z * brho) / (amu * mass)) / \
    np.sqrt(((c * z * brho) / (amu * mass))**2 + 1)

# initial velocity
v0 = np.array([beta * np.sin(a_init / 1000),
               beta * np.cos(a_init / 1000),
               0])


# *******   define magnetic field   *******
def magnetic_field(r):
    x, y, z = r
    b_x = 0
    b_y = 0
    b_z = Bz.BforXplane(x, y)
    return np.array([b_x, b_y, b_z])


# initial conditions of calculation
y0 = np.concatenate([r0, v0])

# *******   motion equation   *******
def lorentz_force(t, y):
    r = y[:3]  # position
    v = y[3:]  # velocity
    B = magnetic_field(r)
    drdt = v
    dvdt = (z*c / (mass*amu*gamma)) * np.cross(v, B)
    return np.concatenate([drdt, dvdt])


# time range of Runge-Kutta
t_span = (0, stop_time/(1e9/c))  # conver time unit in ns

print(f'gamma is {gamma}')
print(f'beta is {beta}')

# *******   solve motion equation using RK45   *******
# conver time unit in ns
solution = solve_ivp(lorentz_force, t_span, y0,
                     method='RK45', max_step=step_time/(1e9/c))

# results
t = solution.t  # timing information
x, y, z = solution.y[0], solution.y[1], solution.y[2]  # positions
vx, vy, vz = solution.y[3], solution.y[4], solution.y[5]  # velocities

# correct revolution time
if y[-1]<0:
    print(f'!!!!!   Break   !!!!!')
    print(f'stop_time ({stop_time}) is too short')
    exit()
else:
    print('******* Revolution time *******')
    print(f'{(t[-1]*1e9/c)-(y[-1]/(vy[-1]/(1e9/c)))} ns')  # conver time unit in ns

print('*******   The number of iterations   *******')
print(f'{len(t)} times')

print('*******   Final positions   *******')
print(f'x is {x[-1]*1000} mm')
print(f'y is {y[-1]*1000} mm')

print('*******   Final velocities   *******')
print(f'vx is {vx[-1]*1000/(1e9/c)} mm/ns')  # conver unit in mm/ns
print(f'vy is {vy[-1]*1000/(1e9/c)} mm/ns')  # conver unit in mm/ns

# *******   plot   *******
# t vs x, and t vs y
box1 = plt.figure(figsize=(15, 5))
fig1_1 = box1.add_subplot(1, 1, 1)
fig1_1 = plt.plot(t/(1e9/c), x*1e3, label='x-motion')
fig1_1 = plt.plot(t/(1e9/c), y*1e3, label='y-motion')
plt.ylabel('position [mm]')
plt.xlabel('time [ns]')
plt.legend()
# x vs y
box2 = plt.figure(figsize=(7, 6))
fig2_1 = box2.add_subplot(1, 1, 1)
fig2_1 = plt.plot(x*1e3, y*1e3, label="X-Y plane")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.axis('equal')
plt.legend()

# stop of an execution timer
end = time.perf_counter()
print('******* Execution time *******')
print(f'{end - start:.2f} seconds')

plt.show()
# plt.savefig('a.pdf')
