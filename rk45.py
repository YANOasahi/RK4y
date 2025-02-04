import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
from fractions import Fraction
from scipy.constants import c
from scipy.constants import physical_constants

import mag_field2 as mag
# import variables as vs
# import position as pos
# import map

# start of an execution timer
start = time.perf_counter()

# *******   physics constant   *******
# light speed [m/s] is defined by scipy.constants
# amu [ev/c^2]
amu = physical_constants['atomic mass unit-electron volt relationship'][0]

# *********************************************
# **************   modify here   **************
# *********************************************
# target nuclide
z, mass = 32.0, 77.922853  # 78Ge, mass is taken from AME2020
# brho of particles when dp/p=0 [Tm]
brho0 = 4.7447
# beta of the ring in the X-axis
betax = 7.817
# the end time of Runge-Kutta
stop_time = 390  # in ns
# stop_time = 378.937  # in ns
# stop_time = 9.7  # in ns
# step time of Runge-Kutta
step_time = 0.1  # max. 100 ps step
# step_time = 0.0001  # max. 100 fs step
# *********************************************

# *******   positions of particles   *******
x0 = 9287.959673
y0 = 0.0
z0 = 15.0
r0 = np.array([x0/1000.0, y0/1000.0, z0/1000.0])  # initial position
print(f'initial position is ({x0}, {y0}, {z0})')

# *******   velocities of particles   *******
# for searching an ideal condition
j, l, k = 0.0, 0.0, 0.0
dp = 0.1 * j
brho = brho0 * (1 + dp / 100.0)
emi = l * 20.0
a_max = np.sqrt(emi / betax) * 0.95
a_bin = 2.0 * a_max / 15.0
a_init = k * a_bin - a_max

# q/m
MoQ = Fraction(mass/z)
# momentum based on brho
mom = Fraction(c*z*brho*1E-6)
# beam energy based on momentum
energy = Fraction((np.sqrt(np.square(mom)+np.square((mass*amu)))-(mass*amu))/mass)
# beta based on brho and q/m
beta = Fraction(((c*z*brho)/(amu*mass))/np.sqrt(((c*z*brho) /
                                        (amu*mass))*((c*z*brho)/(amu*mass))+1.0))
# Lorentz factor
gamma = Fraction(1/np.sqrt(float(1-np.square(beta))))

# initial velocity
v0 = np.array([beta * np.sin(a_init / 1000.0),
               beta * np.cos(a_init / 1000.0),
               0.0])
print(f'initial velocity is ({v0[0]}, {v0[1]}, {v0[2]})')

print('******* Initial conditions *******')
print(f'gamma is {float(gamma):.5f}')
print(f'beta is {float(beta):.5f}')
print(f'beam energy is {float(energy*1e6):.3f} MeV/u')

# notrim_map, trim_map = map.map_making()

# def mag_field(r):
#     # find the nearest magnet
#     x, y, z = r
#     distance_mag=[[((vs.magnet_pos_x[i][j]-(x*1000))**2 + 
#                     (vs.magnet_pos_y[i][j]-(y*1000))**2 +
#                     (vs.magnet_pos_z[i][j]-z*1000)**2) for j in range(4)] for i in range(6)]
#     near_mag = np.unravel_index(np.argmin(distance_mag), (6, 4))
#     min_distance_mag=distance_mag[near_mag[0]][near_mag[1]]

#     # if the nearest magnet is too far,
#     # magnetic field will be np.array([0, 0, 0])
#     if min_distance_mag > 1000**2 + 200**2 + 18**2:
#         print('too far!')
#         return np.array([0, 0, 0])
    
#     elif near_mag[1]==0 or near_mag[1]==3:
#         # convert particle position to the nearest magnet's coordinate
#         particle_pos = pos.Position(x*1000, y*1000, z,
#                                vs.magnet_pos_x[near_mag[0]][near_mag[1]],
#                                vs.magnet_pos_y[near_mag[0]][near_mag[1]],
#                                vs.bend_angle[near_mag[0]][near_mag[1]])
#         # find the nearest map data
#         pos_in_map = [((trim_map[i][0]-particle_pos[0][1])**2 + 
#                        (trim_map[i][1]-particle_pos[0][0])**2 +
#                        (trim_map[i][2]-particle_pos[0][2])**2) 
#                        for i in range(16800000)]
#         nearest = np.argmin(pos_in_map)
#         print('trim magnet')
#         return np.array(
#             [notrim_map[nearest][4], 
#              notrim_map[nearest][3], 
#              notrim_map[nearest][5]]
#             )
    
#     elif near_mag[1]==1 or near_mag[1]==2:
#         # convert particle position to the nearest magnet's coordinate
#         particle_pos = pos.Position(x*1000, y*1000, z,
#                                vs.magnet_pos_x[near_mag[0]][near_mag[1]],
#                                vs.magnet_pos_y[near_mag[0]][near_mag[1]],
#                                vs.bend_angle[near_mag[0]][near_mag[1]])
#         # find the nearest map data
#         pos_in_map = [((notrim_map[i][0]-particle_pos[0][1])**2 + 
#                        (notrim_map[i][1]-particle_pos[0][0])**2 +
#                        (notrim_map[i][2]-particle_pos[0][2])**2) 
#                        for i in range(16800000)]
#         nearest = np.argmin(pos_in_map)
#         print('no trim magnet')
#         return np.array(
#             [trim_map[nearest][4], 
#              trim_map[nearest][3], 
#              trim_map[nearest][5]]
#             )

# initial conditions of calculation
init = np.concatenate([r0, v0])

# *******   motion equation   *******
def lorentz_force(t, y):
    r = y[:3]  # position
    v = y[3:]  # velocity
    B = mag.mag_field(r)
    # differential of position
    drdt = v
    # differential of velocity
    dvdt = Fraction((z*c / (mass*amu*gamma))) * np.cross(v, B)
    return np.concatenate([drdt, dvdt])


# time range of Runge-Kutta
t_span = (0, stop_time/(1e9/c))  # conver time unit in ns

# *******   solve motion equation using RK45   *******
solution = solve_ivp(lorentz_force, t_span, init,
                     method='RK45', max_step=step_time/(1e9/c))

# results
t = solution.t  # timing information
x, y, z = solution.y[0], solution.y[1], solution.y[2]  # positions
vx, vy, vz = solution.y[3], solution.y[4], solution.y[5]  # velocities

# correct revolution time
# if y[-1]<0:
#     print(f'!!!!!   Break   !!!!!')
#     print(f'stop_time ({stop_time}) is too short')
#     exit()
# else:
#     print('******* Revolution time *******')
#     print(f'{(t[-1]*1e9/c)-(y[-1]/(vy[-1]/(1e9/c))):.3f} ns')  # time unit is ns

print('******* Revolution time *******')
print(f'{(t[-1]*1e9/c)-(y[-1]/(vy[-1]/(1e9/c))):.3f} ns')  # time unit is ns

print('*******   The number of iterations   *******')
print(f'{len(t)} times')

print('*******   Final positions   *******')
print(f'x is {x[-1]*1e3} mm')  # convert unit in mm
print(f'y is {y[-1]*1e3} mm')  # convert unit in mm
print(f'z is {z[-1]*1e3} mm')  # convert unit in mm

print('*******   Final velocities   *******')
print(f'vx is {vx[-1]*1000/(1e9/c):.3f} mm/ns')  # convert unit in mm/ns
print(f'vy is {vy[-1]*1000/(1e9/c):.3f} mm/ns')  # convert unit in mm/ns
print(f'vz is {vz[-1]*1000/(1e9/c):.3f} mm/ns')  # convert unit in mm/ns

# *******   output file   *******
with open("rk45_output.dat", "w") as file:
    for column1, column2 in zip(x*1e3, y*1e3):
        file.write(f"{column1} {column2}\n")

# *******   plot   *******
# t vs x, and t vs y
box1 = plt.figure(figsize=(15, 5))
fig1_1 = box1.add_subplot(1, 1, 1)
fig1_1 = plt.plot(t/(1e8/c), x*1e3, label='x-motion')
fig1_1 = plt.plot(t/(1e8/c), y*1e3, label='y-motion')
plt.ylabel('position [mm]')
plt.xlabel('time [ns]')
plt.legend()
# t vs z
box12 = plt.figure(figsize=(15, 5))
fig1_2 = box12.add_subplot(1, 1, 1)
fig1_2 = plt.plot(t/(1e8/c), z*1e3, label='z-motion')
plt.ylabel('position in Z [mm]')
plt.xlabel('time [ns]')
plt.legend()
# x vs y
box2 = plt.figure(figsize=(10.5, 8.4))
fig2_1 = box2.add_subplot(1, 1, 1)
fig2_1 = plt.plot(x*1e3, y*1e3, label="X-Y plane")
abesan=np.genfromtxt('./kidou_long_10ps.dat')
fig2_2 = plt.plot(abesan[:,1], abesan[:,2], label="Abe-san results")
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
